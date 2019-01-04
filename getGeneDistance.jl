#!/usr/bin/env julia
#=
@name: getGeneDistance.jl
@author: Juan C. Castro <jccastrog at gatech dot edu>
@update: 05-Oct-2018
@version: 1.0
@license: GNU General Public License v3.0.
please type "./getGeneDistance.jl -h" for usage help
=#

###===== 1.0 Load packages, define functions, initialize variables =====###
# 1.1 Load packages ======================================================#
using ArgParse;
using GZip;
using NamedTuples;
# 1.2 Define functions ===================================================#
"""
	parse_commandline()
Parse arguments in the command line and output help menu if required 
arguments are not provided.

#Examples

```
julia ./getGeneDistance.jl
required option --ogs_list was not provided
usage: getGeneDistance -l OGS_LIST -g GFF_LIST -f OGS_FILE [-o OUTPUT]

julia ./getGeneDistance.jl -h
usage: getGeneDistance.jl -l OGS_LIST -g GFF_LIST -f OGS_FILE [-o OUTPUT] [-h]

optional arguments:
  -l, --ogs_list OGS_LIST
                        File with the OGs for which distances are
                        going to be calculated.
  -g, --gff_list GFF_LIST
                        File with a list for the location of the gff
                        or gff.gz              files of the genomes
                        included.
  -f, --ogs_file OGS_FILE
                        File with the OGs for the genomes in the
                        pangenome.
  -o, --output OUTPUT   Output file with pairwise distances for the
                        OGs in the list in the             genomes
                        sampled. (default: "ogs_distances.tsv")
  -h, --help            show this help message and exit

"""
function parseCommandline()
    s = ArgParseSettings()

    @add_arg_table s begin
    	"--ogs_list", "-l"
    		help = "File with the OGs for which distances are going to be calculated."
    		required = true
        "--gff_list", "-g"
            help = "File with a list for the location of the gff or gff.gz 
            files of the genomes included."
			required = true
        "--ogs_file", "-f"
            help = "File with the OGs for the genomes in the pangenome."
			required = true
		"--output", "-o"
            help = "Output file with pairwise distances for the OGs in the list in the
            genomes sampled."
			required = false
			default = "ogs_distances.tsv"
    end
    return parse_args(s)
end
                                             """
	getGenomeLength(genome_file::String)
Calculate the length of a given genome if there are multiple contigs 
the length is calculated as the sum of all of the contig lengths.
The return value is and Int64 that denotes the genome length in bp.

#Examples

```julia
makeBlastDB("fasta_file.fa");
28090197

```
"""
function getGenomeLength(genome_file::String)
	seq = "";
	open(genome_file) do fasta_file
		for line in eachline(fasta_file)
			line = rstrip(line);
			if startswith(line, ">")
				continue;
			else
				seq = "$seq$line";
			end
		end
	end
	seq_len = length(seq);
	return(seq_len);
end
"""
	makeBlastDB(fasta_file::String)
Use "makeblastdb" to create nucleotide database files for a fasta file.

#Examples

```julia
makeBlastDB("fasta_file.fa");


Building a new:q DB, current time: 10/05/2018 10:42:49
New DB name:   /home/user/fasta_file
New DB title:  fasta_file.fa
Sequence type: Nucleotide
	Keep MBits: T
Maximum file size: 1000000000B
Adding sequences from FASTA; added 12454 sequences in 0.449869 seconds.
```
"""
function makeBlastDB(fasta_file::String)
	extension(url::String) = try 
		match(r"\.[A-Za-z0-9]+$", url).match 
	catch 
		"" 
	end ;
	fasta_ext = extension(fasta_file);
	basename_file = basename(fasta_file);
	out_file = replace(basename_file, fasta_ext, "");
	blast_cmd = `makeblastdb -in $fasta_file -input_type fasta -dbtype nucl -out $out_file`;
	run(blast_cmd);
end
"""
	executeBlastTbl(query_file::String, db_file::String, num_thread::Int)
Use "blastn" to create map a query to a database and output a table in blast format 6.

#Examples

```julia
executeBlastTbl("query_file.fa", "db_file", 4);
Running blast on 4 threads... Done.

```
"""
function executeBlastTbl(query_file::String, db_file::String, num_thread::Int=3)
	extension(url::String) = try match(r"\.[A-Za-z0-9]+$", url).match catch "" end ;
	query_ext = extension(query_file);
	basename_file = basename(query_file);
	out_file = replace(basename_file, query_ext, "");
	blast_cmd = `blastn -db $db_file -query $query_file -num_threads $num_thread -out $out_file -outfmt 6`;
	print(STDOUT, "Running blast on $num_thread threads... ")
	run(blast_cmd);
	print(STDOUT, "Done.\n")
end
"""
	parseGFF(ggf_file::String, gene_arr::Array)
With a gff or gff.gz file extract the coordinates for each gene as a tuple element
with names start and end.

#Parameters

gff_file File in gff or gff.gz format indicating the coordinates of genes in a given genome.
gene_arr Array with the gene names to get distances for.

#Return

dist_arr Array containing all the pair distances in gene_arr. 

"""
function parseGFF(gff_file::String, gene_arr::Array)
	gene_dict = Dict();
	genome_length = Int64;
	dist_arr = [];
	extension(url::String) = try match(r"\.[A-Za-z0-9]+$", url).match catch "" end ;
	gff_ext = extension(gff_file)
	if gff_ext == ".gz"
		GZip.open(gff_file) do g_file
			for line in eachline(g_file)
				line =  rstrip(line);
				fields = split(line, '\t');
				cds_start = fields[4];
				cds_end = fields[5];
				cds_attr = fields[6];
				attr_fields = split(fields, ";");
				gene_id = split(attr_fields[1], "=")[2];
				if gene_id in gene_arr
					gene_dict["gene_id"] = @NT(s = cds_start, e = cds_end); #s start e end
				end
			end
		end
	elseif gff_ext == ".gff"
		open(gff_file) do g_file
			for line in eachline(g_file)
				line =  rstrip(line);
				fields = split(line, '\t');
				cds_start = fields[4];
				cds_end = fields[5];
				cds_attr = fields[6];
				attr_fields = split(fields, ";");
				gene_id = split(attr_fields[1], "=")[2];
				if gene_id in gene_arr
					gene_dict["gene_id"] = @NT(s = cds_start, e = cds_end); #s start e end
				end
			end
		end
	else
		print(STDERR, "ERROR! File $gff_file does not have a .gff or .gff.gz extension");
	end
	for i in gene_arr
		d1 = Int64;
		d2 = Int64;
		for j in gene_arr
			g1 = i;
			g2 = j;
			g1_start = gene_dict[g1].s;
			g1_end = gene_dict[g1].e;
			g2_start = gene_dict[g2].s;
			g2_end = gene_dict[g2].e;
			if g1_start > g2_start
				d1 = g1_start - g2_end;
				d2 = genome_length - g1_end + g2_start;
			else
				d1 = g2_start - g1_end;
				d2 = genome_length - g1_end + g2_start;
			end
			if d1 < d2
				push!(dist_arr, d1);
			else
				push!(dist_arr, d2);
			end
		end
	end
	return(dist_arr);
end
"""
	parseOGS(ogs_file::String)
With an OGS table extract the sets of genes in each genomes of interest.

#Parameters

ogs_file File where the first row indicates the genome names, and every row
	after indicates an orthologus group (OG). And every colunm represents a genome.

ogs_arr  Array with the OGs to be extracted.

#Return

dist_arr Array containing all the pair distances in gene_arr. 

"""
function parseOGS(ogs_file::String, ogs_arr::Array)

end

parsed_args = parseCommandline();