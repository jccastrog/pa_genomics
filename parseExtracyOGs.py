#!/usr/bin/env python
'''
@name: parseExtracyOGs.py
@author: Juan C. Castro <jccastrog at gatech dot edu>
@update: 05-Jan-2018
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
				og = 'OG_{}'.format(k)
				og_dict[og] = dict()
				for i in xrange(0,num_genomes):
					genome = genome_ids[i]
					if fields[i] == '-':
						continue
					else :
						gene = fields[i].split(',')[0]
						og_dict[og][genome] = gene
				k += 1
	return(k, og_dict)
def  extractOGs(og_dict, fasta_dir):
	for og in og_dict:
		fasta_name = '{}.fasta'.format(og)
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
	num_ogs,ogs_dict = getOGdict(args.ogs_file)
	extractOGs(ogs_dict, args.genome_dir)
'''============================================================'''
