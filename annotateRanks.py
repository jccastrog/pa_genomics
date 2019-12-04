#!/usr/bin/env python
'''
@name: annotateRanks.py
@author: Juan C. Castro <jccastrog at gatech dot edu>
@update: 27-Apr-2018
@version: 1.0
@license: GNU General Public License v3.0.
please type "./annotateRanks.py -h" for usage help
'''

'''===== 1.0 Import modules, define functions, and initialize variables ====='''
#============================= 1.1 Import modules =============================#
try :
	import os, sys
	import argparse
except:
	sys.stderr.write('ERROR! Cannot import required modules (requires os, sys, and argparse)\n')
	sys.exit()
#============================ 1.2 Define functions ============================#
def argParser():
	parser = argparse.ArgumentParser(description=":parseGenebank2tsv: Take a genebank file and extract the function annotations and locus IDs [jccastrog@gatech.edu]")
	group = parser.add_argument_group('Required arguments')
	group.add_argument('-a', action='store', dest='annotations', required=True, help='A file with SEED annotations.')
	group.add_argument('-r', action='store', dest='ranks', required=True, help='A file relating simplified gene names and ortologous groups.')
	args = parser.parse_args()
	return(args)
def annotateRankedOGs(annotations, ranks):
	sys.stdout.write("Rank\tOG\tFunction\n")
	annot_dict = dict()
	with open(annotations) as annotation_file:
		lines = annotation_file.readlines()
		for line in lines:
			if line.startswith("Genome"):
				continue
			else :
				line = line.rstrip("\n")
				fields = line.split("\t")
				OG = fields[4]
				function = fields[5]
				annot_dict[OG] = function
	with open(ranks) as rank_file:
		lines = rank_file.readlines()
		for line in lines:
			if line.startswith("rank"):
				continue
			else :
				line = line.rstrip("\n")
				fields = line.split("\t")
				rank = fields[1]
				OG = fields[2]
				try:
					function = annot_dict[OG]
				except:
					function = "NA"
				out_str = "{}\t{}\t{}".format(rank, OG, function)
				sys.stdout.write(out_str)
				sys.stdout.write("\n")
'''========================== 2.0 Parse the files ==========================='''
if __name__=="__main__":
	args = argParser()
	annotations = args.annotations
	ranks = args.ranks
	annotateRankedOGs(annotations, ranks)
'''=============================================================================='''
