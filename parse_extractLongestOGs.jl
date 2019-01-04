#!/usr/bin/env julia
#=
@name: parse_extractLongestOGs.jl
@author: Juan C. Castro <jccastrog at gatech dot edu>
@update: 25-Dec-2018
@version: 1.0.2
@license: GNU General Public License v3.0.
please type "./parse_extractLongestOGs.jl -h" for usage help
=#
###===== 1.0 Load packages, define functions, initialize variables =====###
# 1.1 Load packages ======================================================#
using ArgParse
# 1.2 Define functions ===================================================#
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--ogsFile", "-f"
            help = "File with the OGs for the genomes in the pangenome."
			required = true
        "--genomeDir", "-d"
            help = "Directory containing the fasta files for the corresponding genomes"
			required = true
		"--extension", "-e"
            help = "The extension of the files in the genomes directory either \"fna\" or \"faa\"."
			required = true
			default = "faa"
		"--output", "-o"
            help = "Output file with sequences for the orthologous groups."
			required = true
			default = "ogs.fasta"
    end
    return parse_args(s)
end
function getDictOG(ogs_file)
	dict_OG = Dict();
	genome_arr = []
	open(ogs_file) do o_file
		k = 0;
		for line in eachline(o_file);
			line = rstrip(line);
			i = 1;
			if k == 0
				genome_arr = split(line, "\t");
				k += 1;
			else
				og_name = "OG_$k";
				dict_OG[og_name] = Dict();
				fields = split(line, "\t");
				for field in fields
					genome_name = genome_arr[i];
					if field == "-"
						i += 1;
						continue;
					else
						genes = split(field, ",");
						dict_OG[og_name][genome_name] = genes;
					end
					i += 1;
				end
				k += 1;
			end
		end
	end
	return genome_arr, dict_OG;
end	
function parseFasta(fasta_file)
	fasta_parse = Dict()
	open(fasta_file) do f_file
		id = "";
		seq = "";
		for line in eachline(f_file)
			line = rstrip(line);
			if startswith(line, ">")
				id = replace(line, ">", "")
				fasta_parse[id] = ""
				seq = ""
			else
				seq = "$seq$line";
				fasta_parse[id] = seq
			end
		end
	end
	return fasta_parse
end
function extractGenes(genome_arr, dict_OG, genome_dir, ext)
	pan_dict = Dict()
	seqs_dict = Dict();
	num_ogs = length(collect(keys(dict_OG)));
	for genome_name in genome_arr
		if ext == "fna"
			fasta_file = "$genome_dir/$genome_name.fna";
		elseif ext == "faa"
			fasta_file = "$genome_dir/$genome_name.faa";
		end
		try
			seqs_dict[genome_name] = parseFasta(fasta_file);
		catch
			write(SDTERR, "ERROR! Cannot find file $genome_name.$ext in directory $genome_dir.\n")
			write(SDTERR, "\tVerify the file exists in the specified directory.\n")
	end
	for ogs in keys(dict_OG)
		seq = "";
		pan_id = "";
		gen_dict = dict_OG[ogs];
		for loc_genome in keys(gen_dict)
			try
				for id in gen_dict[loc_genome]
					loc_seq = seqs_dict[loc_genome][id];
					loc_pan_id = ">$ogs\|$loc_genome\|$id";
					if length(loc_seq) > length(seq)
						seq = loc_seq;
						pan_id = loc_pan_id;
					end
				end
			catch
				print("")
			end

		end
		pan_dict[pan_id] = seq;
		lenSeq = length(seq);
	end
	return pan_dict;
end
# 1.3 Initialize variables ================================================#
parsedArgs = parse_commandline();
ogs_file = parsedArgs["ogsFile"];
genome_dir = parsedArgs["genomeDir"];
output = parsedArgs["output"];
ext = parsedArgs["extension"];
###====================== 2.0 Parse the pangenome ======================###
# 2.1 Parse ogs file =====================================================#
genome_arr, dict_OG = getDictOG(ogs_file);
# 2.2 Extract the genes ==================================================#
pan_dict = extractGenes(genome_arr, dict_OG, genome_dir, ext);
###====================== 2.0 Parse the pangenome ======================###
###======================= 3.0 Write the output ========================###
open(output, "w") do out
	for id in keys(pan_dict)
		seq = pan_dict[id];
		write(out, "$id\n");
		write(out, "$seq\n");
	end
end
#=========================================================================#