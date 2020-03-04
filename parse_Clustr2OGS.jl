#!/usr/bin/env julia
#=
@name: parse_Clustr2OGS.jl
@author: Juan C. Castro <jccastrog at gatech dot edu>
@update: 19-dec-2019
@version: 1.0.2
@license: GNU General Public License v3.0.
please type "./parse_Clustr2OGS.jl -h" for usage help
=#

###===== 1.0 Load packages, define functions, initialize variables =====###
# 1.1 Load packages ======================================================#
using ArgParse;
# 1.2 Define functions ===================================================#
function parseCommandline()
	s = ArgParseSettings();
	
	@add_arg_table s begin
		"--clstr_file", "-c"
			help = "Output file with the clustering of CD-HIT";
			required = true;
		"--strain_file", "-s"
			help = "File with the names of the strains used";
			required = true;
		"--output", "-o"
			help = "File name for the output";
			required = false;
			default = "output.tsv";
	end
return parse_args(s);
end
function parseClsrtFile(clstr::String)
	k = 0;
	ogs_dict = Dict();
	strain_dict = Dict();
	open(clstr) do file
		for line in eachline(file)
			if startswith(line, ">clustering")
				k += 1;
				strain_dict = Dict()
			else
				fields = split(line, " ");
				og = split(fields[2], "|");
				strain = replace(og[1], ">" => "");
				gene = replace(og[2], "..." => "");
				if haskey(strain_dict, strain)
					push!(strain_dict[strain], gene);
				else
					strain_dict[strain] = [gene];
				end
				ogs_dict["OG_$k"] = strain_dict;
			end
		end
	end
	return(ogs_dict)
end
function parseStrainFile(strain_file::String)
	strain_list = [];
	open(strain_file) do sf
		for line in eachline(sf)
			line = rstrip(line);
			push!(strain_list, line);
		end
	end
	return(strain_list);
end
function writeTable(ogs_dict::Dict, strain_list::Array,output::String)
	out_file = open(output, "w");
	for strain in strain_list
		write(out_file, "\t$strain");
	end
	write(out_file, "\n");
	for key in collect(keys(ogs_dict))
		write(out_file, "$key\t");
		for strain in strain_list
			if haskey(ogs_dict[key], strain)
				if length(ogs_dict[key][strain]) > 1
					write(out_file, join(ogs_dict[key][strain], ","))
				else
					write(out_file, "$(ogs_dict[key][strain][1])");
				end
				write(out_file, "\t");
			else
				write(out_file, "-\t");
			end
		end
		write(out_file, "\n");
	end
	close(out_file);
end
parsedArgs = parseCommandline();
clstr_file = parsedArgs["clstr_file"];
strain_file = parsedArgs["strain_file"];
output = parsedArgs["output"];
ogs_dict = parseClsrtFile(clstr_file);
strain_list = parseStrainFile(strain_file);
writeTable(ogs_dict, strain_list,  output);
###########################################################################