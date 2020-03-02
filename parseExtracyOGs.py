#!/usr/bin/env python
'''
@name: parseExtracyOGs.py
@author: Juan C. Castro <jccastrog at gatech dot edu>
@update: 02-Mar-2020
@version: 1.0
@license: GNU General Public License v3.0.
please type "./parseExtracyOGs.py -h" for usage help
'''

'''===== 1.0 Import modules, define functions, and initialize variables ====='''
#===== 1.1 Import modules =====
try :
	import os, sys
	import argparse
	from Bio import SeqIO
except:
	sys.stderr.write('ERROR! Cannot import required modules (requires os, sys, argparse, time, progressbar, and Bio)\n')
	sys.exit()

#===== 1.2 Define functions =====
def argParser():
	parser = argparse.ArgumentParser(description="parseExtracyOGs: ID [jccastrog@gatech.edu]")
	group = parser.add_argument_group('Required arguments')
	group.add_argument('-f', action='store', dest='ogs_file', required=True, help='A file with the ogs for the genomes in the pangneome.')
	group.add_argument('-d', action='store', dest='genome_dir', required=True, help='A directory including the fasta files for the corresponding genomes.')
	group.add_argument('-o', action='store', dest='out_dir', required=True, help='Directory where the output file will be stored.')
	args = parser.parse_args()
	return(args)
def getOGdict(ogs_file):
	og_dict = dict()
	with open(ogs_file) as ogs:
		k = 0
		lines = ogs.readlines()
		for line in lines:
			line = line.rstrip('\n')
			if k == 0:
				header = line
				genome_ids = header.split('\t')
				num_genomes = len(genome_ids)
				k += 1
			else :
				fields = line.split('\t')
				#og = 'OG_{}'.format(k) # For miga formating
				og = fields[0] # For cd-hit clustr formating
				og_dict[og] = dict()
				for i in xrange(1,num_genomes):
					genome = genome_ids[i]
					if fields[i] != '-':
						gene = fields[i].split(',')[0]
						#print(genome+"\t"+gene)
						og_dict[og][genome] = gene
						#print(og_dict[og])
				k += 1
	return(k, og_dict)
def  extractOGs(og_dict, fasta_dir, out_dir):
	for og in og_dict:
		fasta_name = '{}/{}.fasta'.format(out_dir,og)
		out_file = open(fasta_name, 'wb')
		current_dict = og_dict[og]
		for genome in current_dict:
			genome_file = '{}/{}.fna'.format(fasta_dir,genome)
			#genome_file = '{}/{}.faa'.format(fasta_dir, genome)
			with open(genome_file) as fasta_file:
				for fasta_parse in SeqIO.parse(fasta_file,"fasta"):
					ID = fasta_parse.id
					seq = fasta_parse.seq
					if current_dict[genome] == ID:
						id_name = '>{}|{}\n'.format(genome, ID)
						out_file.write(id_name)
						out_file.write(str(seq)+'\n')
		out_file.close()

'''================== 2.0 Download the files =================='''
if __name__=="__main__":
	args = argParser()
	num_ogs,ogs_dict = getOGdict(args.ogs_file)
	extractOGs(ogs_dict, args.genome_dir, args.out_dir)
'''============================================================'''
