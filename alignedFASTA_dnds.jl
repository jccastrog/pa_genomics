#!/usr/bin/env julia
#=
@name: alignedFASTA_dnds.jl
@author: Juan C. Castro <jccastrog at gatech dot edu>
@update: 02-Mar-2020
@version: 1.0.0
@license: GNU General Public License v3.0.
please type "./alignedFASTA_dnds.jl -h" for usage help
=#

###===== 1.0 Load packages, define functions, initialize variables =====###
# 1.1 Load packages ======================================================#
using ArgParse;
# 1.2 Define functions ===================================================#
"""
	parse_commandline()
Parse arguments in the command line and output help menu if required 
arguments are not provided.

#Examples

```
julia ./alignedFASTA_dnds.jl
required option --pangenome_file was not provided
usage: getGeneDistance -p PANGENOME_FILE -g PANGENOME_FILE -r REF_FILE -c CONFIG_FILE [-o OUTPUT]

julia ./alignedFASTA_dnds.jl -h
usage: alignedFASTA_dnds.jl -p PANGENOME_FILE -g CDS_DIR -r REF_FILE
                      -c CONFIG_FILE [-n NUM_THREADS] [-o OUTPUT] [-h]

optional arguments:
  -a, --alignment_file	ALIGNMENT_FILE
  						File with aligned sequences for a given orthologous
                        group.
  -o, --output OUTPUT   Output file. Three columns are printed first columns
                        is DN values second column is DS values final column
                        is DN/DS.
  -h, --help            show this help message and exit
```
"""
function parseCommandline()
	s = ArgParseSettings();
	
	@add_arg_table s begin
	"--alignment_file", "-c"
	help = "File with aligned sequences for a given orthologous
            group.";
	required = true;
	"--output", "-o"
	help = "Output file. Three columns are printed first columns
            is DN values second column is DS values final column
            is DN/DS.";
	required = false;
	default = "output.tsv";
	end
	return parse_args(s);
end
"""
	splitByCodon(seq::String)
Split a nucleotide sequence for a protein in triplets representing codons.

#Input
	seq::Sring => A nucleotide sequence for a protein.
	
#Return
	codons::Array => An array where each element is a codon sequence.

#Examples
```
julia> seq1 = splitByCodon("ATGAAACCCGGGTTTTAA")
6-element Array{Any,1}:
 "ATG"
 "AAA"
 "CCC"
 "GGG"
 "TTT"
 "TAA"

```
"""
function splitByCodon(seq::String)
	codons = [];
	l_seq = length(seq);
	if l_seq%3==0
		for i in 1:Int(l_seq/3)
				codon = seq[i*3-2:i*3];
				push!(codons, codon);
		end
	end
	return codons;
end
"""
	singleDnDs(seq1::Array, seq2::Array, genetic_code::Dict, n_dict::Dict, s_dict::Dict)
Calculate Dn/Ds values for a pair of orthologous aligned seuqneces using 
	the Nei-Gojobori method (Nei & Gojobori, 1986).

#Input
	seq1::Array => An array with a nucleotide sequence split by codons.
	
	seq2::Array => An array with a nucleotide sequence split by codons.
	
	genetic_code => A dictionary with the genetic code translating codons to
					aminoacids.

	n_dict => A dictionary with the expected values for nonsynonymous subsitutions.

	s_dict => A dictionary with the expected values for synonymous subsitutions.

#Return
	dn_ds::Tuple => A tuple with values for DN and DS of the two aequences being compared.

#Examples
```
julia> singleDnDs(seq1, seq2, genetic_code, n_dict, s_dict)
(dn = 0.1505, ds = 1.2074)

```
"""
function singleDnDs(seq1::Array, seq2::Array, genetic_code::Dict, n_dict::Dict, s_dict::Dict)
	N = 0;
	S = 0;
	Nd = 0
	Sd = 0;
	l_seq1 = length(seq1);
	l_seq2 = length(seq2);
	if l_seq1==l_seq2
		for i in 1:l_seq1
			N += n_dict[seq1[i]];
			S += s_dict[seq1[i]];
			if seq1[i]!= seq2[i]	
				if genetic_code[seq1[i]] != genetic_code[seq2[i]];
					Nd += 1;
				else
					Sd += 1;
				end

			end
		end
		pn = Nd/N;
		ps = Sd/S;
		dn = (-3/4)*log(1- (4*pn/3));
		ds = (-3/4)*log(1- (4*ps/3));
		dn_ds = (dn = dn, ds = ds);
	end
	return(dn_ds);
end
"""
	parseFasta(fasta_file::String)
In a fasta file of aligned sequences extract the sequences into a dictionary.

#Input
	fasta_file::String => Path to the fasta file with aligned sequences.
	
#Return
	seq_dict::Dict => Dictionary matching sequence identifiers and sequences.

"""
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
"""
	allDnDs(fasta_file::String, genetic_code::Dict, n_dict::Dict, s_dict::Dict)
In a fasta file of aligned sequences calculate the DN and DS values for each
	pair of sequences.

#Input
	fasta_file::String => Path to the fasta file with aligned sequences.
	
	genetic_code => A dictionary with the genetic code translating codons to
					aminoacids.

	n_dict => A dictionary with the expected values for nonsynonymous subsitutions.

	s_dict => A dictionary with the expected values for synonymous subsitutions.

#Return
	arr_dn_ds::Array => Array with values for DN and DS.

"""
function allDnDs(fasta_file::String, genetic_code::Dict, n_dict::Dict, s_dict::Dict)
	seq_dict = parseFasta(fasta_file);
	seq_ids = collect(keys(seq_dict))
	arr_dnds = [];
	i = 0;
	j = 0;
	for s1 in seq_ids
		i += 1;
		if !occursin("-", s1)
			seq1 = splitByCodon(s1);
			for s2 in seq_ids
				j += 1;
				if !occursin("-", s2)
					seq2 = splitByCodon(s2)
					if j > i
						loc_dnds = singleDnDs(seq1, seq2, genetic_code, n_dict, s_dict);
						push!(arr_dnds, loc_dnds);
					end
				end
			end
		end
	end
	return(arr_dnds);
end
# 1.3 Initialize variables ===============================================#
# 1.3.1 Codon variables #
genetic_code = Dict("TTT" => "Phe",
			 		"TTC" => "Phe",
			 		"TTA" => "Leu",
			 		"TTG" => "Leu",
			 		"CTT" => "Leu",
			 		"CTC" => "Leu",
			 		"CTA" => "Leu",
			 		"CTG" => "Leu",
			 		"ATT" => "Ile",
			 		"ATC" => "Ile",
			 		"ATA" => "Ile",
			 		"ATG" => "Met",
			 		"GTT" => "Val",
			 		"GTC" => "Val",
			 		"GTA" => "Val",
			 		"GTG" => "Val",
			 		"TCT" => "Ser",
			 		"TCC" => "Ser",
			 		"TCA" => "Ser",
			 		"TCG" => "Ser",
			 		"CCT" => "Pro",
			 		"CCC" => "Pro",
			 		"CCA" => "Pro",
			 		"CCG" => "Pro",
			 		"ACT" => "Thr",
			 		"ACC" => "Thr",
			 		"ACA" => "Thr",
			 		"ACG" => "Thr",
			 		"GCT" => "Ala",
			 		"GCC" => "Ala",
			 		"GCA" => "Ala",
			 		"GCG" => "Ala",
			 		"TAT" => "Tyr",
			 		"TAC" => "Tyr",
			 		"TAA" => "Stop",
			 		"TAG" => "Stop",
			 		"CAT" => "His",
			 		"CAC" => "His",
			 		"CAA" => "Gln",
			 		"CAG" => "Gln",
			 		"AAT" => "Asn",
			 		"AAC" => "Asn",
			 		"AAA" => "Lys",
			 		"AAG" => "Lys",
			 		"GAT" => "Asp",
			 		"GAC" => "Asp",
			 		"GAA" => "Glu",
			 		"GAG" => "Glu",
			 		"TGT" => "Cys",
			 		"TGC" => "Cys",
			 		"TGA" => "Stop",
			 		"TGG" => "Trp",
			 		"CGT" => "Arg",
			 		"CGC" => "Arg",
			 		"CGA" => "Arg",
			 		"CGG" => "Arg",
			 		"AGT" => "Ser",
			 		"AGC" => "Ser",
			 		"AGA" => "Arg",
			 		"AGG" => "Arg",
			 		"GGT" => "Gly",
			  		"GGC" => "Gly",
			  		"GGA" => "Gly",
			  		"GGG" => "Gly"
);
n_dict = Dict("TTT" => 2+2/3,
			  "TTC" => 2+2/3,
			  "TTA" => 2,
			  "TTG" => 2,
			  "CTT" => 2,
			  "CTC" => 2,
			  "CTA" => 1+1/3,
			  "CTG" => 1+1/3,
			  "ATT" => 2+1/3,
			  "ATC" => 2+1/3,
			  "ATA" => 2+1/3,
			  "ATG" => 3,
			  "GTT" => 2,
			  "GTC" => 2,
			  "GTA" => 2,
			  "GTG" => 2,
			  "TCT" => 2,
			  "TCC" => 2,
			  "TCA" => 2,
			  "TCG" => 2,
			  "CCT" => 2,
			  "CCC" => 2,
			  "CCA" => 2,
			  "CCG" => 2,
			  "ACT" => 2,
			  "ACC" => 2,
			  "ACA" => 2,
			  "ACG" => 2,
			  "GCT" => 2,
			  "GCC" => 2,
			  "GCA" => 2,
			  "GCG" => 2,
			  "TAT" => 2+2/3,
			  "TAC" => 2+2/3,
			  "TAA" => 2+1/3,
			  "TAG" => 2+2/3,
			  "CAT" => 2+2/3,
			  "CAC" => 2+2/3,
			  "CAA" => 2+2/3,
			  "CAG" => 2+2/3,
			  "AAT" => 2+2/3,
			  "AAC" => 2+2/3,
			  "AAA" => 2+2/3,
			  "AAG" => 2+2/3,
			  "GAT" => 2+2/3,
			  "GAC" => 2+2/3,
			  "GAA" => 2+2/3,
			  "GAG" => 2+2/3,
			  "TGT" => 2+2/3,
			  "TGC" => 2+2/3,
			  "TGA" => 2+2/3,
			  "TGG" => 3,
			  "CGT" => 2,
			  "CGC" => 2,
			  "CGA" => 1+2/3,
			  "CGG" => 1+2/3,
			  "AGT" => 2+2/3,
			  "AGC" => 2+2/3,
			  "AGA" => 2+1/3,
			  "AGG" => 2+1/3,
			  "GGT" => 2,
			  "GGC" => 2,
			  "GGA" => 2,
			  "GGG" => 2
);
s_dict = Dict("TTT" => 1/3,
			  "TTC" => 1/3,
			  "TTA" => 1,
			  "TTG" => 1,
			  "CTT" => 1,
			  "CTC" => 1,
			  "CTA" => 1+2/3,
			  "CTG" => 1+2/3,
			  "ATT" => 2/3,
			  "ATC" => 2/3,
			  "ATA" => 2/3,
			  "ATG" => 0,
			  "GTT" => 1,
			  "GTC" => 1,
			  "GTA" => 1,
			  "GTG" => 1,
			  "TCT" => 1,
			  "TCC" => 1,
			  "TCA" => 1,
			  "TCG" => 1,
			  "CCT" => 1,
			  "CCC" => 1,
			  "CCA" => 1,
			  "CCG" => 1,
			  "ACT" => 1,
			  "ACC" => 1,
			  "ACA" => 1,
			  "ACG" => 1,
			  "GCT" => 1,
			  "GCC" => 1,
			  "GCA" => 1,
			  "GCG" => 1,
			  "TAT" => 1/3,
			  "TAC" => 1/3,
			  "TAA" => 2/3,
			  "TAG" => 1/3,
			  "CAT" => 1/3,
			  "CAC" => 1/3,
			  "CAA" => 1/3,
			  "CAG" => 1/3,
			  "AAT" => 1/3,
			  "AAC" => 1/3,
			  "AAA" => 1/3,
			  "AAG" => 1/3,
			  "GAT" => 1/3,
			  "GAC" => 1/3,
			  "GAA" => 1/3,
			  "GAG" => 1/3,
			  "TGT" => 1/3,
			  "TGC" => 1/3,
			  "TGA" => 1/3,
			  "TGG" => 0,
			  "CGT" => 1,
			  "CGC" => 1,
			  "CGA" => 1+1/3,
			  "CGG" => 1+1/3,
			  "AGT" => 1/3,
			  "AGC" => 1/3,
			  "AGA" => 2/3,
			  "AGG" => 2/3,
			  "GGT" => 1,
			  "GGC" => 1,
			  "GGA" => 1,
			  "GGG" => 1
);
# 1.2.2 Parser variables #
parsed_args = parseCommandline();
alignment_file = parsed_args["alignment_file"];
output = parsed_args["output"];
###======== 2.0 Calculate DN/DS for protein sequences in a file ========###
arr_dnds = allDnDs(alignment_file, genetic_code, n_dict, s_dict);
# 1.3 Initialize variables ===============================================#
# 2.1 Parse the array and print the values to STDOUT =====================#
for i in arr_dnds
	loc_tup = arr_dnds[i];
	loc_dnds = loc_tup.dn/loc_tup.ds;
	wrt_line = "$(loc_tup.dn)\t$(loc_tup.ds)\t$loc_dnds\n";
	write(STDOUT, wrt_line);
end
###=====================================================================###