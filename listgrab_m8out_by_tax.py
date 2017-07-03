#!/usr/bin/env python

from os import sys

# usage: listgrabpy listfile csvfile
SeqListFile = sys.argv[1]
seqlist = open(SeqListFile, 'r')
inpath = sys.argv[2]
infile = open(inpath, 'r')
outfile_path = sys.argv[3]

# load outfile
OutFile = open(outfile_path,'w')

getset=set([])	# declare list
for line in seqlist:
	line = line.strip()
	getset.add(line)	# append line to list
# print getset

# example:
# S23C1_TRINITY_DN1607639_c7_g1_i1_2_1,Florenciella_sp_Strain_RCC1587_tax236786_locID_CAMPEP_0182558112_seqID14928359,236786,370.2,-99.2

# we want to: split the file and keep the contig_id [0] and the tax_id [2] (actually we just need taxid! we can split contig_id with awk -F, later
# if the tax_id is in the getset, then output the line
for line in infile:
	tax_id = line.split(',')[2]
	if tax_id in getset:
		OutFile.write(line)

OutFile.close()
