#!/usr/bin/env julia
#=
@name: parse_concatenateGenes.jl
@author: Juan C. Castro <jccastrog at gatech dot edu>
@update: 04-Sep-2019
@version: 1.0.0
@license: GNU General Public License v3.0.
please type "./parse_concatenateGenes.jl -h" for usage help
=#

###===== 1.0 Load packages, define functions, initialize variables =====###
# 1.1 Load packages ======================================================#
using ArgParse;
# 1.2 Define functions ===================================================#
function parseCommandline()
	s = ArgParseSettings();
	
	@add_arg_table s begin
	"--fasta_directory", "-f"
	help = "The directory storing the .fasta files";
	required = true;
	"--ogs_ids", "-i"
	help = "The file relating gene names with orthlogous groups";
	required = true;
	"--output", "-o"
	help = "File name for the output";
	required = false;
end
return parse_args(s);
end
function parseFasta(fasta_file::String)
	seq_dict = Dict();seq = "";
	id = "";
	sed_line = true;
	open(fasta_file) do f
		for line in eachline(f)
			line = rstrip(line);
			if startswith(line, ">");
				fields = split(line, "|");
				id = fields[1];
				id = replace(id, ">" => "");
				seq = "";
			else
				seq = "$seq$line";
				seq_dict["$id"] = seq;
			end
		end
	end
	return(seq_dict);
end
function parseIDs(ids_file::String)
	ess_genes = [];
	open(ids_file) do ids
		for line in eachline(ids)
			line = rstrip(line);
			fields = split(line, "\t");
			push!(ess_genes, fields[1]);
		end
	end
	return(ess_genes);
end
function concatGenes(fasta_directory::String, ids_file::String, output::String)
	ess_genes = parseIDs(ids_file);
	concat_seqs = Dict();
	for (root, dirs, files) in walkdir(fasta_directory)
		for file in files
			if endswith(file, "fasta")
				og_id = split(file, ".")[1];
				if og_id in ess_genes
					seq_dic = parseFasta(joinpath(root,file));
					for key in collect(keys(seq_dic))
						if haskey(concat_seqs, key)
							old_seq = concat_seqs[key];
							new_seq = seq_dic[key];
							concat_seqs[key] = "$old_seq$new_seq";
						else
							concat_seqs[key] = seq_dic[key];
						end
					end
				end
			end
		end
	end
	out_file = open(output, "w");
	for key in collect(keys(concat_seqs))
		write(out_file, ">$key\n");
		write(out_file, "$(concat_seqs[key])\n");
	end
	close(out_file);
end
function main()
	parsedArgs = parseCommandline();
	fasta_directory = parsedArgs["fasta_directory"];
	ids_file = parsedArgs["ogs_ids"];
	output = parsedArgs["output"];
	concatGenes(fasta_directory, ids_file, output);
end
main();
