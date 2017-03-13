#!/usr/bin/env python

# dedupe_fasta_by_defline.py

"""Reads in a FASTA file. Will deduplicate by DEFLINE, not sequence content.
Basically, this will output (keep) the first incidence of a sequence with a unique defline
and discard any further matches. If this is output of hmm, the order is in descending quality
so the best seq with any given defline will be kept.
If input is gene.fasta, outputs to gene.dd.fasta.
"""
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("input_fasta", help="FASTA-formatted input file")
args = parser.parse_args()

def build_output_handle(infile_path):
	handle_elts = infile_path.split(".")
	handle_elts.insert(-1,"dd")
	outfile_path = ".".join(handle_elts)
	return outfile_path

# load the input fasta and output fasta
input_fasta = open(args.input_fasta, 'r')
outfile_path = build_output_handle(args.input_fasta)
output_fasta = open(outfile_path, 'w')

DeflineDict = {}

for seq_record in SeqIO.parse(input_fasta, "fasta"):
	seq_id = seq_record.id.split("/")[0]
	if seq_id not in DeflineDict:
		DeflineDict[seq_id] = 1
		SeqIO.write(seq_record, output_fasta, "fasta")

input_fasta.close()
output_fasta.close()
