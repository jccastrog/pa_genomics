#!/usr/bin/env julia
#=
@name: getGeneAges.jl
@author: Juan C. Castro <jccastrog at gatech dot edu>
@update: 23-Jul-2019
@version: 1.1
@license: GNU General Public License v3.0.
please type "./getGeneAges.jl -h" for usage help
=#

###===== 1.0 Load packages, define functions, initialize variables =====###
# 1.1 Load packages ======================================================#
using ArgParse;
using Distributed;
# 1.2 Initialize variables ===============================================#
# 1.2.1 Create a parser #
"""
	parse_commandline()
Parse arguments in the command line and output help menu if required 
arguments are not provided.

#Examples

```
julia ./getGeneAges.jl
required option --pangenome_file was not provided
usage: getGeneDistance -p PANGENOME_FILE -g PANGENOME_FILE -r REF_FILE -c CONFIG_FILE [-o OUTPUT]

julia ./getGeneAges.jl -h
usage: getGeneAges.jl -p PANGENOME_FILE -g CDS_DIR -r REF_FILE
                      -c CONFIG_FILE [-n NUM_THREADS] [-o OUTPUT] [-h]

optional arguments:
  -p, --pangenome_file PANGENOME_FILE
                        Pangenome file in fasta format where each
                        sequence is an orthologous group.
  -g, --cds_dir CDS_DIR
                        Directory where the CDS files for each
                        genome are stored in fasta format.
  -r, --ref_file REF_FILE
                        File with a list in which each entry is a
                        filename and the corresponding taxnonomy.
  -c, --config_file CONFIG_FILE
                        Configuration file, this file contains the
                        relevant taxonomic levels used to identify
                        the gene ages.
  -n, --num_threads NUM_THREADS
                        Number of threads to use. (default: 3)
  -o, --output OUTPUT   Output file where the age of each gene is
                        displayed. (default: "output.tsv")
  -h, --help            show this help message and exit
```
"""
function parseCommandline()
	s = ArgParseSettings();
	
	@add_arg_table s begin
		"--pangenome_file", "-p"
			help = "Pangenome file in fasta format where each
			sequence is an orthologous group.";
			required = true;
		"--cds_dir", "-g"
			help = "Directory where the CDS files for each
			genome are stored in fasta format.";
			required = true;
		"--ref_file", "-r"
			help = "File with a list in which each entry is a
			filename and the corresponding taxnonomy.";
			required = true;
		"--config_file", "-c"
			help = "Configuration file, this file contains the
			relevant taxonomic levels used to identify
			the gene ages.";
			required = true;   	
		"--num_threads", "-n"
			help = "Number of threads to use."
			required = false;
			default = 3;
		"--output", "-o"
			help = "Output file where the age of each gene is
			displayed.";
			required = false;
			default = "output.tsv";
	end
	return parse_args(s);
end
# 1.2.2 Parser variables #
parsed_args = parseCommandline();
pangenome_file = parsed_args["pangenome_file"];
cds_dir = parsed_args["cds_dir"];
ref_file = parsed_args["ref_file"];
config_file = parsed_args["config_file"];
num_threads = parse(Int, parsed_args["num_threads"]);
output = parsed_args["output"];
# 1.2.3 Global variables #
addprocs(num_threads);
num_threads_blast = 3;
min_iden = float(40);
# 1.3 Define functions ===================================================#
"""
	makeBlastDB(fasta_file::String)
Use "makeblastdb" to create nucleotide database files for a fasta file.

#Input
	fasta_file::String => A file in fasta format from which a database will
	be built.
#Return
	out_file::String => The path to the output of the database, for future
		   reference.

#Examples

```julia


Building a new DB, current time: 10/05/2018 10:42:49
New DB name:   /home/user/fasta_file
New DB title:  fasta_file.fa
Sequence type: Nucleotide
	Keep MBits: T
Maximum file size: 1000000000B
Adding sequences from FASTA; added 12454 sequences in 0.449869 seconds.

/home/user/fasta_file
```
"""
function makeBlastDB(fasta_file::String, out_dir::String)
	extension(url::String) = try
		match(r"\.[A-Za-z0-9]+$", url).match 
	catch 
		"" 
	end ;
	fasta_ext = extension(fasta_file);
	basename_file = basename(fasta_file);
	out_file = replace(basename_file, fasta_ext => "");
	out_file = "$out_dir/$out_file";
	blast_cmd = `makeblastdb -in $fasta_file -input_type fasta -dbtype prot -out $out_file`;
	run(blast_cmd);
	return(out_file);
end
"""
         parseConfigFile(config_file::String)
Create a directory that specifies the name of each taxonomic rank based on 
the configuration file.

#Inputs
	config_file::String => File defining the names for each taxonomic rank.
#Return
	config_dict::Dict => A dictionary that relates taxonomic level and name.
#Examples
```julia
parseConfigFile("config_file.txt")

Dict{String,String} with 7 entries:
  "phylum"  => "Proteobacteria"
  "family"  => "Enterobacteriaceae"
  "Species" => "Escherichia_coli"
  "genus"   => "Escherichia"
  "domain"  => "Bacteria"
  "class"   => "Gammaproteobacteria"
  "order"   => "Enterobacteriales"
```
For a config file with the format
d: Bacteria
p: Proteobacteria
c: Gammaproteobacteria
o: Enterobacteriales
f: Enterobacteriaceae
g: Escherichia
s: Escherichia_coli

"""
function parseConfigFile(config_file::String)
	config_dict = Dict();
        open(config_file, "r") do cf
                for line in eachline(cf)
                        line = rstrip(line);
                        fields = split(line, " ");
                        if fields[1] == "d:"
                                config_dict["domain"] = fields[2];
                        elseif fields[1] == "p:"
                                config_dict["phylum"] = fields[2];
                        elseif fields[1] == "c:"
                                config_dict["class"] = fields[2];
                        elseif fields[1] == "o:"
                                config_dict["order"] = fields[2];
                        elseif fields[1] == "f:"
                                config_dict["family"] = fields[2];
                        elseif fields[1] == "g:"
                                config_dict["genus"] = fields[2];
                        elseif fields[1] == "s:"
                                config_dict["species"] = fields[2];
                        end
                end
        end
	return config_dict;
end
"""
         createStrataDBs(config_file::String, ref_file::String, cds_dir)
Create a blast database for each stratum (taxonomic rank). Separating the files in the
cds directory into ranks and then concatenating them. Finally use makeBlastDB to create 
separate databases.

#Inputs
	config_file::String => File defining the names for each taxonomic rank.
	cds_dir::String => Directory wehre the CDS of each genome are stored in fasta 
:q
aminoacid format.
	ref_file::String => File defining the detailed taxonomy for each of the CDS files.
#Return
	database_paths::Dict => An dictionary where each element denotes the path of a
			temporary database for each taxonomic rank.
"""
function createStrataDBs(config_file::String, ref_file::String, cds_dir::String)
#This could be optimized with pmap
	database_paths = Dict();
	config_dict = parseConfigFile(config_file);
	domain_name = config_dict["domain"];
	phylum_name = config_dict["phylum"];
	class_name = config_dict["class"];
	order_name = config_dict["order"];
	family_name = config_dict["family"];
	genus_name = config_dict["genus"];
	species_name = config_dict["species"];
	fasta_dir = mktempdir();
	parseConcatFasta(source::String, destination::String) = open(source) do s
		d = open(destination, "a");
		bn = basename(source);
		for line in eachline(s)
			if startswith(line, ">")
				write(d, "$line|$bn\n");
			else
				write(d, "$line\n");
			end
		end
		close(d);
	end
	open(ref_file) do rf
		for line in eachline(rf)
			line = rstrip(line);
			fields = split(line, "\t");
			taxonomy = split(fields[2], " ");
			filename = fields[1];
			filename = "$filename.faa";
			filename = "$cds_dir/$filename";
			if length(taxonomy) >= 7
				if isfile(filename)
					if taxonomy[7] == "s:$species_name"
						parseConcatFasta(filename, "$fasta_dir/species.faa");
					elseif taxonomy[6] == "g:$genus_name"
						parseConcatFasta(filename, "$fasta_dir/genus.faa");
					elseif taxonomy[5] == "f:$family_name"
						parseConcatFasta(filename, "$fasta_dir/family.faa");
					elseif taxonomy[4] == "o:$order_name"
						parseConcatFasta(filename, "$fasta_dir/order.faa");
					elseif taxonomy[3] == "c:$class_name"
						parseConcatFasta(filename, "$fasta_dir/class.faa");
					elseif taxonomy[2] == "p:$phylum_name"
						parseConcatFasta(filename, "$fasta_dir/phylum.faa");
					elseif taxonomy[1] == "d:$domain_name"
						parseConcatFasta(filename, "$fasta_dir/domain.faa");
					end
				end
			end
		end
	end
	for (root, dirs, files) in walkdir(fasta_dir)
		for file in files
			db_name = "$fasta_dir/$file";
			db_path = makeBlastDB(db_name, fasta_dir);
			tax_rank = basename(db_path);
			database_paths[tax_rank] = db_path;
		end
	end
	return database_paths;
end
@everywhere begin
	"""
		executeBlastTbl(query_file::String, db_file::String, num_thread::Int)
	Use "blastn" to create map a query to a database and output a table in blast format 6.
	
	#Inputs
		query_file::String => A file path of the query.
		db_file::String => A file path for the database.
		num_thread::Int => The number of threads to use (default: 3).
	#Return
		out_file::String => Path to the output table for future reference
	Examples
	
	``julia
	executeBlastTbl("query_file.fa", "db_file", 4);
	Running blast for db_file on 4 threads... Done.
	
	db_file.tbl
	```
	"""
	function executeBlastTbl(blast_args::Dict)
		query_file = blast_args["query"];
		db_file = blast_args["db"];
		num_thread = try blast_args["num_thread"]
		catch
			3
		end
		extension(url::String) = try match(r"\.[A-Za-z0-9]+$", url).match 
		catch
			""
		end;
		basename_db = basename(db_file);
		out_file = "$basename_db.tbl";
		blast_cmd = `blastp -db $db_file -query $query_file -num_threads $num_thread -out $out_file -outfmt 6 -num_alignments 1`;
		write(stdout, "Running BLAST for $basename_db on $num_thread threads... ");
		run(blast_cmd);
		write(stdout, "Done.\n");
		return out_file;
	end
end
"""
	parallelBlastTbl(query_file::String, db_dict::Dict, num_thread_blast::Int, num_threads::Int)
Use executeBlastTbl to map sequneces in a query to multiple databases in parallel, each mapping is done
separately and has its own multi thread settings.

#Inputs
	query_file: A file path to the query, a single query is used for all of the blasts.
	db_dict: An array where the key is a taxonomic rank and the value is the path
		 for a blast database.
	num_thread_blast: Number of processors to use for a single blast.
	num_threads: Number of htreads on which to paralellize the blast submissions.
#Return
	tbl_paths: A dictionary with the paths for each of the output tables resulting of the blast,
		   where the key is a taxonomic rank and the value is the file path.
"""
function parallelBlastTbl(query_file::String, db_dict::Dict, num_thread_blast::Int=3)
	input_args = [];
	tax_ranks = collect(keys(db_dict));
	for tr in tax_ranks
		db_file = db_dict[tr];
		blast_args = Dict("query" => query_file, "db" => db_file, "num_thread" => num_thread_blast);
		push!(input_args, blast_args);
	end
	num_workers = nprocs() - 1;
	write(stdout, "Running parallel BLAST in $num_workers workers\n");
	tbl_paths = pmap(executeBlastTbl, input_args);
	write(stdout, "BLAST finished for all ranks!\n");
	return(tbl_paths);
end
"""
	filterBlastTables(tbl_paths::Dict, min_aln::Int, output::String)
Take the files resulting from parallelBlastTbl and orderly assing strata by filtering
the basal strata first moving towards the more derived taxonomic rank.

#Inputs
	tbl_paths:  A dictionary with the paths for each of the output tables resulting of the blast,
                    where the key is a taxonomic rank and the value is the file path.
	min_aln: Minimum percentage of alignment for agene to be considered orthologous,
			is calculated based on the length of the query sequence.
	min_iden_length: Minimum percentage of alignment identity to be considered an ortholog.
"""
function filterBlastTables(tbl_paths::Array, min_iden::Float64, output::String)
	num_tables = length(tbl_paths);
	write(stdout, "Parsing $num_tables to assign taxonomic ranks...");
	assigned_seqs = Dict();
	parse_tbl_file(tbl_path::String, assigned_seqs::Dict, tax_rank::String) = open(tbl_path) do df
		for line in eachline(df)
			fields = split(line, "\t");
			query_id = fields[1];
			subject_id = fields[2]
			pident = parse(Float64,fields[3]);
			if pident >= min_iden
				if haskey(assigned_seqs, query_id)
					""
				else
					assigned_seqs[query_id] = tax_rank;
				end
			end
		end
	end
	try
		tbl_path = tbl_paths[Bool[occursin("domain",i) for i in tbl_paths]][1]
		parse_tbl_file(tbl_path, assigned_seqs, "domain");
	catch
		write(stderr, "Warning!, no sequences mapped to domain stratum!\n");
	end
	try
		tbl_path = tbl_paths[Bool[occursin("phylum",i) for i in tbl_paths]][1]
		parse_tbl_file(tbl_path, assigned_seqs, "phylum");
	catch
		write(stderr, "Warning!, no sequences mapped to phylum stratum!\n");
	end
	try
		tbl_path = tbl_paths[Bool[occursin("class",i) for i in tbl_paths]][1]
		parse_tbl_file(tbl_path, assigned_seqs, "class");
	catch
		write(stderr, "Warning!, no sequences mapped to class stratum!\n");
	end
	try
		tbl_path = tbl_paths[Bool[occursin("order",i) for i in tbl_paths]][1]
		parse_tbl_file(tbl_path, assigned_seqs, "order");
	catch
		write(stderr, "Warning!, no sequences mapped to order stratum!\n");
	end
	try
		tbl_path = tbl_paths[Bool[occursin("family",i) for i in tbl_paths]][1]
		parse_tbl_file(tbl_path, assigned_seqs, "family");
	catch
		write(stderr, "Warning!, no sequences mapped to family stratum!\n");
	end
	try
		tbl_path = tbl_paths[Bool[occursin("genus",i) for i in tbl_paths]][1]
		parse_tbl_file(tbl_path, assigned_seqs, "genus");
	catch
		write(stderr, "Warning!, no sequences mapped to genus stratum!\n");
	end
	try
		tbl_path = tbl_paths[Bool[occursin("species",i) for i in tbl_paths]][1]
		parse_tbl_file(tbl_path, assigned_seqs, "species");
	catch
		write(stderr, "Warning!, no sequences mapped to species stratum!\n");
	end
	out_file = open(output, "w");
	for (seq, stratum) in assigned_seqs
		write(out_file, "$seq\t$stratum\n");
	end
	close(out_file);
	write(stdout,"Done!\nTaxonomic rank assignments have been saved to $output\n");
ends
###========== 2.0  Create databases and perform the mappings ===========###
# 2.1 Create databases ===================================================#
database_paths = createStrataDBs(config_file, ref_file, cds_dir);
# 2.2 Map OGs recursively ================================================#
tbl_paths = parallelBlastTbl(pangenome_file, database_paths, num_threads_blast);
###============ 3.0  Filter mapping and assign gene strata =============###
filterBlastTables(tbl_paths, min_iden, output);
###=====================================================================###
