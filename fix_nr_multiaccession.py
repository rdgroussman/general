#!/usr/bin/env python

# fix messy genbank headers
# Ryan Groussman
# Armbrust lab
# September 2016

import sys

filename = sys.argv[1]

fasta_file = open(filename,'r')
OutFileName = filename + ".filtered"
OutFile = open(OutFileName,'w')

for line in fasta_file:
    if line.startswith(">"):
        line_minus = line[1:]
        cutoff = line_minus.find(">")
        outline = line[:cutoff] + "\n"
        OutFile.write(outline)
        print "a header..."
    else:
        OutFile.write(line)
        print "not a header..."
OutFile.close()
