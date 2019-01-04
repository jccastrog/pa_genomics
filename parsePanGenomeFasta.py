#!/usr/bin/env python
'''
@name: parsePanGenomeFasta.py
@author: Juan C. Castro <jccastrog at gatech dot edu>
@update: 05-Jan-2017
@version: 1.0
@license: GNU General Public License v3.0.
please type "./parsePanGenomeFasta.py -h" for usage help
'''

'''===== 1.0 Import modules, define functions, and initialize variables ====='''
#===== 1.1 Import modules =====
try :
	import os, sys
	import argparse
	import time, progressbar
	from Bio import SeqIO
except:
	sys.stderr.write('ERROR! Cannot import required modules (requires os, sys, argparse, time, progressbar, and Bio)\n')
	sys.exit()

#===== 1.2 Define functions =====
def argParser():
	parser = argparse.ArgumentParser(description=":parsePanGenomeFasta: Take a representative for each OG in the pangenome and add it to a fasta file [jccastrog@gatech.edu]")
	group = parser.add_argument_group('Required arguments')
	group.add_argument('-f', action='store', dest='ogs_file', required=True, help='A file with the ogs for the genomes in the pangneome.')
	group.add_argument('-d', action='store', dest='genome_dir', required=True, help='A directory including the fasta files for the corresponding gneomes.')
	group.add_argument('-o', action='store', dest='out_file', default = 'ogs.fasta', required=True, help='Output file with sequences for the orthologous groups. (default : %(default)s)')
	args = parser.parse_args()
	return(args)
def getGeneList(ogs_file):
	gene_dict = dict()
	with open(ogs_file) as ogs:
		lines = ogs.readlines()
		k = 0
		for line in lines:
			line = line.rstrip('\n')
			if k == 0:
				genome_list = line.split('\t')
				k += 1
			else :
				fields = line.split('\t')
				i = 0
				for field in fields:
					if len(field.split(',')) > 1:
						i += 1
						continue
					else:
						if field == '-':
							i += 1
							continue
						else:
							if genome_list[i] in gene_dict:
								gene_dict[genome_list[i]].append(field)
							else:
								gene_dict[genome_list[i]] = [field]
							i += 1
							break
				k += 1
	return(gene_dict)	
def createPanFasta(fasta_dir, gene_dict, fasta_out):
	num_genes = len(gene_dict)
	out_file = open(fasta_out, 'wb')
	with progressbar.ProgressBar(max_value=num_genes) as bar:
		i = 1
		for file in os.listdir(fasta_dir):
			if file.endswith('.fna') or file.endswith('fa'):
				genome_name = os.path.splitext(file)[0]
				if genome_name in gene_dict:
					with open(file) as fasta_file:
						print len(SeqIO.parse(fasta_file, "fasta"))
						for fasta_parse in SeqIO.parse(fasta_file, "fasta"):
							ID = fasta_parse.id
							seq = fasta_parse.seq
							if ID in gene_dict[genome_name]:
								out_file.write('>{}|{}\n'.format(genome_name,ID))
								out_file.write('{}\n'.format(seq))
								bar.update(i)
								i += 1
	out_file.close()

'''================== 2.0 Download the files =================='''
if __name__=="__main__":
	args = argParser()
	gene_dict = getGeneList(args.ogs_file)
	createPanFasta(args.genome_dir, gene_dict, args.out_file)
'''============================================================'''
