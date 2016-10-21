#!/usr/bin/env python


"""
Extracts sequences from a FASTA file given a list of MMETSP IDs.

Example:
MMETSP1428
MMETSP1429
MMETSP1065
MMETSP0580
"""

from Bio import SeqIO
import argparse

# declare argument parser stuff
parser = argparse.ArgumentParser()
parser.add_argument("input_fasta", help="FASTA-formatted input file")
parser.add_argument("input_list", help="List of MMETSP ids to keep")
parser.add_argument("output_fasta", help="Name of FASTA-formatted output file")
args = parser.parse_args()

id_list = open(args.input_list, 'r')
fasta = open(args.input_fasta, 'r')
outfile = open(args.output_fasta, 'w')

mmetsp_set = set([])
for mmetsp_id in id_list:
    mmetsp_set.add(mmetsp_id.strip())
id_list.close()

for seq in SeqIO.parse(fasta, "fasta"):
    if seq.id[:10] in mmetsp_set:
        SeqIO.write(seq, outfile, "fasta")
