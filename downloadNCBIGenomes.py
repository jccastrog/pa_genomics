#!/usr/bin/env python
'''
@name: downloadNCBIGenomes.py
@author: Juan C. Castro <jccastrog at gatech dot edu>
@update: 04-Aug-2017
@version: 1.0
@license: GNU General Public License v3.0.
please type "./downloadNCBIGenomes.py -h" for usage help
'''

'''===== 1.0 Import modules, define functions, and initialize variables ====='''
#===== 1.1 Import modules =====
try :
	import os, sys, subprocess
	import argparse
	import gzip, wget
except:
	sys.stderr.write('ERROR! Cannot import required modules (requires os, sys, subprocess, argparse, gzip)\n')
	sys.exit()

#===== 1.2 Define functions =====
def arg_parser():
	parser = argparse.ArgumentParser(description="downloadNCBIGenomes: Download NCBI genomes for a particuar taxon ID [jccastrog@gatech.edu]")
	group = parser.add_argument_group('Required arguments')
	group.add_argument('-t', action='store', dest='taxon_id', required=True, help='The taxon ID of which genomes will be downloaded')
	args = parser.parse_args()
	return(args)
def download_summary():
	if not os.path.exists('assembly_summary.txt'):
		wget.download('ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt')
def download_files(genomeFile, taxid):
	with open(genomeFile) as genomeList:
		lines = genomeList.readlines()
		for line in lines:
			if line.startswith('#'):
				continue
			else:
				fields = line.split('\t')
				spTaxID = fields[6]
				ftpName = fields[19].split('/')[-1]
				ftpName = '{0}/{1}_genomic.fna.gz'.format(fields[19],ftpName)
				if taxid==spTaxID:
					fileName =  os.path.basename(ftpName)
					if not os.path.exists(fileName):
						wget.download(ftpName)

'''===== 2.0 Download the files ====='''
if __name__=="__main__":
	args = arg_parser()
	download_summary()
	summary_file = 'assembly_summary.txt'
	download_files(summary_file, args.taxon_id)

''' ===================='''









