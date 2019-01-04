#!/usr/bin/env python
'''
@name: parseOGS.py
@author: Juan C. Castro <jccastrog at gatech dot edu>
@update: 27-Apr-2018
@version: 1.0
@license: GNU General Public License v3.0.
please type "./parseOGS.py -h" for usage help
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
	group.add_argument('-k', action='store', dest='kbase_annot', required=True, help='A file with SEED annotations.')
	group.add_argument('-o', action='store', dest='og_file', required=True, help='A file relating simplified gene names and ortologous groups.')
	group.add_argument('-r', action='store', dest='rel_file', required=True, help='A file that relates genomes genes and KBase IDs.')
	args = parser.parse_args()
	return(args)
def fileParser(kbase_annot, og_file, rel_file):
	sys.stdout.write("Genome\tGene\tMergedName\tKBaseID\tOG\tFunction\n")
	kbase_dict = dict()
	og_dict = dict()
	with open(kbase_annot) as kbase_file:
		lines = kbase_file.readlines()
		for line in lines:
			line = line.rstrip("\n")
			fields = line.split("\t")
			pID = fields[0]
			function = fields[1]
			kbase_dict[pID] = function
	with open(og_file) as o_file:
		lines = o_file.readlines()
		for line in lines:
			line = line.rstrip("\n")
			fields = line.split("\t")
			simpName = fields[0]
			simpName = simpName.replace("|", "-")
			ogID = fields[1]
			og_dict[simpName] = ogID
	with open(rel_file) as r_file:
		lines = r_file.readlines()
		for line in lines:
			line = line.rstrip("\n")
			fields = line.split("\t")
			genome = fields[0]
			gene = fields[1]
			pID = fields[2]
			simpName = "{}-{}".format(genome, gene)
			ogID = og_dict[	simpName]
			function = kbase_dict[pID]
			outStr = "{}\t{}\t{}\t{}\t{}\t{}".format(genome, gene, simpName, pID, ogID, function)
			sys.stdout.write(outStr)
			sys.stdout.write("\n")
'''========================== 2.0 Parse the files ==========================='''
if __name__=="__main__":
	args = argParser()
	kbase_annot = args.kbase_annot
	og_file = args.og_file
	rel_file = args.rel_file
	fileParser(kbase_annot, og_file, rel_file)
'''=============================================================================='''
