#!/usr/bin/env python

from os import sys

# usage: listgrabpy listfile csvfile
SeqListFile = sys.argv[1]
seqlist = open(SeqListFile, 'r')
inpath = sys.argv[2]
infile = open(inpath, 'r')

# load outfile
OutFileName = SeqListFile + ".contigs.csv"
OutFile = open(OutFileName,'w')

getset=set([])	# declare list
for line in seqlist:
	line = line.strip()[6:]
	getset.add(line)	# append line to list
# print getset

for line in infile:
	seq_id = line.split(',')[0]
	if seq_id in getset:
		OutFile.write(line)

OutFile.close()
