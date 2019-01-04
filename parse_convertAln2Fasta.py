#!/usr/bin/env python
'''
@name: parse_convertAln2Fasta.py
@author: Juan C. Castro <jccastrog at gatech dot edu>
@update: 15-Oct-2018
@version: 1.0.4
@license: GNU General Public License v3.0.
please type "parse_convertAln2Fasta.py -h" for usage help
'''

###===== 1.0 Load packages, define functions, initialize variables =====###
# 1.1 Load packages ======================================================#
import sys, os
import argparse
from Bio import SeqIO
from Bio import AlignIO
# 1.2. Define functions ==================================================
def argParser():
	parser = argparse.ArgumentParser(description="parseExtracyOGs: ID [jccastrog@gatech.edu]")
	group = parser.add_argument_group('Required arguments')
	group.add_argument('-a', action='store', dest='aln_file', required=True, help='Alignment file in clustal format (.aln).')
	group.add_argument('-f', action='store', dest='fasta_file', required=True, help='Alignment file in fasta format (.fasta)..')
	args = parser.parse_args()
	return(args)
def convert_aln_phy(inFile, outFile):
	inHndl = open(inFile, 'r')
	outHndl = open(outFile, 'w')
	inSeqs = [x for x in SeqIO.parse (inFile, 'clustal')]
	alignment = AlignIO.read(inHndl, "clustal")
	outHndl.write(alignment.format("fasta"))

###============================= 2.0  Main ============================###
if __name__=="__main__":
	args = argParser()
	convert_aln_phy(args.aln_file, args.fasta_file)
###====================================================================###
