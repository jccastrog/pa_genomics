#!/usr/bin/env python
'''
@name: parseGenebank2tsv.py
@author: Juan C. Castro <jccastrog at gatech dot edu>
@update: 27-Apr-2018
@version: 1.0
@license: GNU General Public License v3.0.
please type "./parseGenebank2tsv.py -h" for usage help
'''

'''===== 1.0 Import modules, define functions, and initialize variables ====='''
#============================= 1.1 Import modules =============================#
try :
	import os, sys
	import argparse
	import time
	from Bio import SeqIO
	from Bio.Seq import Seq
	from Bio.SeqRecord import SeqRecord
	from Bio.SeqFeature import SeqFeature, FeatureLocation
	from Bio import SeqIO
except:
	sys.stderr.write('ERROR! Cannot import required modules (requires os, sys, argparse, time, and Bio)\n')
	sys.exit()

#============================ 1.2 Define functions ============================#
def argParser():
	parser = argparse.ArgumentParser(description=":parseGenebank2tsv: Take a genebank file and extract the function annotations and locus IDs [jccastrog@gatech.edu]")
	group = parser.add_argument_group('Required arguments')
	group.add_argument('-g', action='store', dest='genebank_file', required=True, help='A genebank file that contains the IDs and annotated annotations for pangenome genes.')
	args = parser.parse_args()
	return(args)
def parseGeneBankAnnot(genebank_file):
	sys.stdout.write("pseudomonasID\tfunction\n")
	recs = [rec for rec in SeqIO.parse(genebank_file, "genbank")]
	for rec in recs:
		feats = [feat for feat in rec.features if feat.type == "gene"]
		for feat in feats:
			pID = feat.qualifiers.items()[0][1][0]
			try :
				function = feat.qualifiers.items()[1][1][0]
			except :
				function = "NA"
			fields = "{}\t{}".format(pID,function)
			sys.stdout.write(fields)
			sys.stdout.write("\n")
'''========================== 2.0 Parse the files ==========================='''
if __name__=="__main__":
	args = argParser()
	parseGeneBankAnnot(args.genebank_file)
'''============================================================'''