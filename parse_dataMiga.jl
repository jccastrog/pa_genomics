#!/usr/bin/env julia
#=
@name: parse_dataMiga.jl
@author: Juan C. Castro <jccastrog at gatech dot edu>
@update: 15-Jan-2018
@version: 1.0.4
@license: GNU General Public License v3.0.
please type "./parse_dataMiga.jl -h" for usage help
=#

###===== 1.0 Load packages, define functions, initialize variables =====###
# 1.1 Load packages ======================================================#
using ArgParse
using ConfParser
# 1.2 Define functions ===================================================#
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--log_directory", "-l"
            help = "The directory storing the .log files"
			required = true
        "--output", "-o"
            help = "File name for the output"
			required = true
    end
    return parse_args(s)
end
function parse_logs(logDir)

	for (root, dirs, files) in walkdir(logDir)
		for file in files
			genomeName = split(file, ".")[1];
			numContigs = 0;
			N50 = 0;
			totalLength = 0;
			proteins = 0;
			averageLength = 0;
			completeness = 0
			contamination = 0;
			quality = 0 ;
			lFIle = open(file)
			for line in eachline(lFIle)
				fields = split(line, " ")
				#print("$fields[1]\n")
				if fields[1]=="Contigs:"
					numContigs = fields[2];
				elseif fields[1]=="N50:"
					N50 = fields[2];
				elseif fields[1]=="Total"
					totalLength = fields[3];
				elseif fields[1]=="Predicted"
					proteins = fields[3];
				elseif fields[1]=="Average"
					averageLength = fields[3];
				elseif fields[1]=="Completeness:"
					completeness = fields[2];
				elseif fields[1]=="Contamination:"
					contamination = fields[2];
				elseif fields[1]=="Quality:"
					quality = fields[2];
				end
			end
			print("$genomeName\t$numContigs\t$N50\t$totalLength\t$proteins\t")
			print("$averageLength\t$completeness\t$contamination\t$quality\n")
		end
	end
end

function main()
	parsedArgs = parse_commandline();
	logDirectory = parsedArgs["log_directory"];
	output = parsedArgs["output"];
	parse_logs(logDirectory);
end

main();
#=
for (root, dirs, files) in walkdir(".")
    println("Directories in $root")
    for dir in dirs
        println(joinpath(root, dir)) # path to directories
    end
    println("Files in $root")
    for file in files
        println(joinpath(root, file)) # path to files
    end
end
=#