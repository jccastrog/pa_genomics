#!/usr/bin/env julia
#=
@name: parse_CATGenesAllGenomes.jl
@author: Juan C. Castro <jccastrog at gatech dot edu>
@update: 21-Oct-2020
@version: 1.1
@license: GNU General Public License v3.0.
please type "./parse-
CATGenesAllGenomes.jl -h" for usage help
=#

###===== 1.0 Load packages, define functions, initialize variables =====###
# 1.1 Load packages ======================================================#
using ArgParse;
# 1.2 Initialize variables ===============================================#
# 1.2.1 Create a parser #
"""
EDIT THIS LATER SERIOUSLY
"""
function parseCommandline()
	s = ArgParseSettings();
	
	@add_arg_table s begin
	"--cds_dir", "-c"
	help = "Directory where the CDS fasta files are 
	stored.";              
	required = true;
	"--ext", "-e"
	help = "Extension of the fasta files
	genome are stored in fasta format.";
	required = true;
	default = ".faa"
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
cds_dir = parsed_args["cds_dir"];
ext = parsed_args["ext"];
output = parsed_args["output"];
# 1.3 Define functions ===================================================#
function parseFasta(fasta_file::String)
	extension(url::String) = try
		match(r"\.[A-Za-z0-9]+$", url).match 
	catch 
		"" 
	end ;
	fasta_ext = extension(fasta_file);
	basename_file = basename(fasta_file);
	genome_name = replace(basename_file, fasta_ext => "");
	seq_dict = Dict();
	seq = "";
	id = "";
	sed_line = true;
	open(fasta_file) do f
		for line in eachline(f)
			line = rstrip(line);
			if startswith(line, ">");
				fields = split(line, "|");
				fields[1] = replace(fields[1], ">" => "");
				id = "$genome_name|$(fields[1])";
				seq = "";
			else
				seq = "$seq$line";
				seq_dict["$id"] = seq;
			end
		end
	end
	return(seq_dict);
end
function parseCDS(cds_dir::String, output::String, ext::String)
	out_file = open(output, "w");
	extension(url::String) = try
		match(r"\.[A-Za-z0-9]+$", url).match 
	catch 
		"" 
	end ;
	for (root, dirs, files) in walkdir(cds_dir)
		for file in files
			file_ext = extension(file);
			if file_ext == ext
				seq_dict = parseFasta(joinpath(root,file));
				for key in collect(keys(seq_dict))
					wrt_line = 	">$key\n$(seq_dict[key])\n";
					write(out_file, wrt_line);
				end
			end
		end
	end
end
parseCDS(cds_dir, output, ext);