#!/usr/bin/env python

from os import sys
from Bio import SeqIO

# usage: listgrabpy listfile fastafile
# load master FASTA file, list of seq names
# eg; Thaps genes: /Users/rgroussman/Dropbox/Armbrust/share/data/seq/organisms/Thalassiosira_pseudonana/genome/
SeqListFile = sys.argv[1]
seqlist = open(SeqListFile, 'r')
FASTAfile = sys.argv[2]
fasta = open(FASTAfile, 'r')

# load outfile
OutFileName = SeqListFile + ".fasta"
OutFile = open(OutFileName,'w')

getlist=[]	# declare list
for line in seqlist:
	line = line.strip()
	getlist.append(line)	# append line to list
print getlist

for record in SeqIO.parse(fasta, 'fasta'):
	if record.id in getlist:
		SeqIO.write(record, OutFile, "fasta")

OutFile.close()
