#!/usr/bin/env python
'''
@name: parseConcatGenes.py
@author: Juan C. Castro <jccastrog at gatech dot edu>
@update: 28-May-2018
@version: 1.0
@license: GNU General Public License v3.0.
please type "./parseConcatGenes.py -h" for usage help
'''

'''===== 1.0 Import modules, define functions, and initialize variables ====='''
#============================= 1.1 Import modules =============================#
try :
	import os, sys
	import argparse, progressbar, time
	from Bio import SeqIO
except:
	sys.stderr.write('ERROR! Cannot import required modules (requires os, sys, argparse, time, progressbar, and Bio)\n')
	sys.exit()
#============================ 1.2 Define functions ============================#
def argParser():
	parser = argparse.ArgumentParser(description="parseConcatGenes: In a directory of files with genes run through the file and concatenate the genes in each file and then create a single fasta file where each entry is a list of concatenated genes. [jccastrog@gatech.edu]")
	group = parser.add_argument_group('Required arguments')
	group.add_argument('-d', action='store', dest='dir_fasta', required=True, help='The directory with the fasta files.')
	group.add_argument('-e', action='store', dest='extension', required=True, default='.fasta', help='The extension of the file to be used. (default : %(default)s)')
	group.add_argument('-o', action='store', dest='output', required=True, default='output.fasta', help='The output file, a fasta file where each entry is a set of concatenated genes from a genome. (default : %(default)s)')
	args = parser.parse_args()
	return(args)
def concatFile(gene_file):
	ID = os.path.splitext(gene_file)[0]
	ext = os.path.splitext(gene_file)[-1]
	big_seq = ''
	with open(gene_file) as fasta_file:
		for fasta_parse in SeqIO.parse(fasta_file,"fasta"):
			seq = fasta_parse.seq
			big_seq = big_seq+seq
	return [ID, big_seq]
def createBigFasta(directory, extension, output):
	files = os.listdir(directory)
	bar = progressbar.ProgressBar()
	output_file = open(output, 'w')
	for file in bar(files):
		if file.endswith(extension):
			filename = '{}/{}'.format(directory, file)
			[ID, seq] = concatFile(filename)
			output_file.write('>{}\n'.format(ID))
			output_file.write('{}\n'.format(seq))
		else :
			continue
def main():
	args = argParser()
	createBigFasta(args.dir_fasta, args.extension, args.output)
'''== 2.0 Run through the files in the directory and concatenate the genes =='''
#================================ 2.1 Run main ================================#
if __name__=="__main__":
	main()
'''=========================================================================='''