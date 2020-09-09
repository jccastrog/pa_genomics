#!/usr/bin/env julia
#=
@name: getGeneDistance.jl
@author: Juan C. Castro <jccastrog at gatech dot edu>
@update: 09-Sep-2020
@version: 1.0
@license: GNU General Public License v3.0.
please type "./getGeneDistance.jl -h" for usage help
=#

###===== 1.0 Load packages, define functions, initialize variables =====###
# 1.1 Load packages ======================================================#
using ArgParse;
using GZip;
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
usage: getGeneDistance.jl -l OGS_LIST -g GFF_LIST -f OGS_FILE
                        [-o OUTPUT] [-h]

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

    @add_arg_table! s begin
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
	extension(url::String) = try 
		match(r"\.[A-Za-z0-9]+$", url).match 
	catch 
		""
	end ;
	gff_ext = extension(gff_file)
	if gff_ext == ".gz"
		GZip.open(gff_file) do g_file
			genome_name = splitext(splitext(gff_file)[1])[1];
			for line in eachline(g_file)
				line =  rstrip(line);
				if startswith(line, "#")
					if startswith(line, "# Sequence Data:")
						fields = split(line, " ");
						length_field = split(fields[4],";")[2];
						genome_length = parse(Int64, split(length_field, "=")[2]);
					end
				else
					fields = split(line, '\t');
					cds_start = parse(Int64, fields[4]);
					cds_end = parse(Int64, fields[5]);
					cds_attr = fields[9];
					attr_fields = split(cds_attr, ";");
					gene_id = split(attr_fields[1], "=")[2];
					gene_dict[gene_id] = (s = cds_start, e = cds_end); #s start e end
				end
			end
		end
	elseif gff_ext == ".gff3"
		open(gff_file) do g_file
			genome_name = splitext(gff_file)[1];
			for line in eachline(g_file)
				line =  rstrip(line);
				if startswith(line, "#")
					if startswith(line, "# Sequence")
						fields = split(line, " ");
						length_field = split(fields[4],";")[2];
						genome_length = parse(Int64, split(length_field, "=")[2]);
					end
				else
					fields = split(line, "\t");
					cds_start = parse(Int64, fields[4]);
					cds_end = parse(Int64, fields[5]);
					cds_attr = fields[9];
					attr_fields = split(cds_attr, ";");
					gene_id = split(attr_fields[1], "=")[2];
					gene_dict[gene_id] = (s = cds_start, e = cds_end); #s start e end
				end
			end
		end
	else
		write(stderr, "ERROR! File $gff_file does not have a .gff3 or .gff3.gz extension");
	end
	arr_length = length(gene_arr);
	for i in 1:arr_length
		d1 = Int64;
		d2 = Int64;
		for j in 1:arr_length
			if i>j
				tuple_1 = gene_arr[i];
				tuple_2 = gene_arr[j];
				for gi in tuple_1.genes
					for gj in tuple_2.genes
						gi_start = gene_dict[gi].s;
						gi_end = gene_dict[gi].e;
						gj_start = gene_dict[gj].s;
						gj_end = gene_dict[gj].e;
						if gi_start > gj_start
							d1 = gi_start - gj_end;
							d2 = genome_length - gi_end + gj_start;
						else
							d1 = gj_start - gi_end;
							d2 = genome_length - gj_end + gi_start;
						end
						if d1 < d2
							loc_arr = [tuple_1.og tuple_2.og d1];
							dist_arr = [dist_arr ; loc_arr];
						else
							loc_arr = [tuple_1.og tuple_2.og d2];
							dist_arr = [dist_arr ; loc_arr];
						end
					end
				end
			end
		end
	end
	return(dist_arr);
end
"""
	parseOGS(ogs_file::String, ogs_arr::Array)
With an OGS table extract the sets of genes in each genomes of interest.

#Parameters

ogs_file File where the first row indicates the genome names, and every row
	after indicates an orthologus group (OG). And every colunm represents a genome.

ogs_arr  Array with the OGs to be extracted.

#Return

dist_arr Array containing all the pair distances in gene_arr. 

	"""
function parseOGS(ogs_file::String, ogs_arr::Array)
	genome_dict = Dict();
	genome_arr = [];
	open(ogs_file) do file
		for line in eachline(file)
			line = rstrip(line);
			if startswith(line, "OG")
				fields = split(line, "\t");
				len_fields = length(fields);
				og = fields[1];
				gene_ids = fields[2:len_fields];
				num_genes = length(gene_ids);
				local_dict = Dict();
				if og in ogs_arr
					for i in 1:num_genes
						if gene_ids[i]!="-"
							single_genes = split(gene_ids[i], ",");
							tuple_genes = (og = og, genes = single_genes)
							local_genome = genome_arr[i];
							if haskey(genome_dict,local_genome)
								genome_dict[local_genome] = vcat(genome_dict[local_genome], tuple_genes);
							else
								genome_dict[local_genome] = tuple_genes;
							end
						end
					end
				end
			else
				genome_arr = split(line, "\t");
				len_arr = length(genome_arr);
				genome_arr = genome_arr[2:len_arr];
				for i in genome_arr
					genome_dict[i] = [];
				end
			end
		end
	end
	return(genome_dict);
end
"""
	getGeneDistances(genome_dict::Dict, gff_list::String)
Use a dictionary that relates genomes to genes to calculate distances between genes
i a list of gff3 files.

#Parameters

genome_dict A genome directoty as created by parseOGS.

gff_list  A file where each row is an absolute path to a gff file.

#Return

all_distances An array contning the distances for thegenes in all the genomes.

"""
function getGeneDistances(genome_dict::Dict, gff_list::String)
	gff_arr = [];
	all_distances = [];
	extension(url::String) = try 
		match(r"\.[A-Za-z0-9]+$", url).match 
	catch 
		""
	end ;
	open(gff_list) do gff_file
		for line in eachline(gff_file)
			line = rstrip(line);
			line = convert(String, line)
			bn_gff = basename(line);
			push!(gff_arr, line);
		end
	end
	for gff_file in gff_arr
		bn_gff = basename(gff_file);
		genome_name = "";
		if extension(gff_file) == ".gff3"
			genome_name = replace(bn_gff, ".gff3" => "");
		elseif  extension(gff_file) == ".gz"
			genome_name = replace(bn_gff, ".gff3.gz" => "");
		end
		gene_arr = genome_dict[genome_name];
		local_distances = parseGFF(gff_file, gene_arr);
		all_distances = vcat(all_distances, local_distances)
	end
	return(all_distances)
end
# 1.3 Initialize variables ===============================================#
# 1.3.1 Parser variables #
parsed_args = parseCommandline();
ogs_list = parsed_args["ogs_list"];
gff_list = parsed_args["gff_list"];
ogs_file = parsed_args["ogs_file"];
output = parsed_args["output"];
# 1.3.2 Global variables #
ogs_arr = [];
###=========== 2.0 Calculate gene distances for a set of OGS ===========###
# 2.1 Parse the OGS file into an array ===================================#
open(ogs_list) do o_l
	for line in eachline(o_l)
		line = rstrip(line);
		push!(ogs_arr, line);
	end
end
# 2.2 Parse the OGS file to relate the gene IDs of each genome ===========#
genome_dict = parseOGS(ogs_file, ogs_arr);
# 2.3 Calculate distances ================================================#
all_distances = getGeneDistances(genome_dict, gff_list);
# 2.3 Calculate distances ================================================#
###======================= 3.0 Write the output ========================###
out_stream = open(output, "w");
for distance_value in all_distances
	write(out_stream, "$distance_value\n");
end
close(out_stream);
#=========================================================================#