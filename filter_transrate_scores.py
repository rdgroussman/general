#!/usr/bin/env python

# filter_transrate_scores.py


"""
Transrate scores are output in 'contigs.csv' which looks something like this:
contig_name,length,prop_gc,orf_length,in_bridges,p_good,p_bases_covered,p_seq_true,score,p_not_segmented,eff_length,eff_count,tpm,coverage,sCnuc,sCcov,sCord,sCseg
TRINITY_DN84_c0_g1_i1,363,0.38292,57,0,0,0.0,1,0.01,0.997424,363,0.0,0.0,0.0,1,0.0,0,0.997424
TRINITY_DN6_c0_g1_i1,333,0.51952,109,0,0,0.444444,1,0.01,0.997424,333,0.0,0.0,0.0,1,0.444444,0,0.997424
TRINITY_DN82_c0_g1_i1,320,0.4125,72,0,0,0.0,1,0.01,0.997424,320,0.0,0.0,0.0,1,0.0,0,0.997424
...
(another 2 million lines)

The task of this script is to parse these scores by desired critera (score, sCseg, etc) and output a new contigs.csv
with all contigs passing the passed criteria.

"""

import argparse
import re

# from Bio import Phylo

parser = argparse.ArgumentParser()
parser.add_argument("input_csv", help="Contigs.csv")
parser.add_argument("-o", "--out_file", help="Output file name", type=str)
parser.add_argument("-s", "--score", help="Lowest acceptable contig score (float between 0 and 1)", type=float)
parser.add_argument("-c", "--sCseg", help="Lowest acceptable Cseg score (float between 0 and 1)", type=float)
args = parser.parse_args()

outfile = open(args.out_file, 'w')

print "Opening contigs.csv..."
input_csv = open(args.input_csv, 'r')
header = next(input_csv)
outfile.write(header)

for line in input_csv:
	line_elts = line.split(",")
	contig_name = line_elts[0]
	score = float(line_elts[8])
	sCseg = float(line_elts[17])
	if score >= args.score and sCseg >= args.sCseg:
		outfile.write(line)

input_csv.close()
outfile.close()
