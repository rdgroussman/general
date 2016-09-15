#!/usr/bin/env python

"""
Ryan Groussman
Armbrust Lab
September 2016

This script takes as input a fasta file and a length argument. Will output a fasta file with all sequences
equal or larger than the specified length.
"""

# import packages
import argparse
from Bio import SeqIO

def build_output_handle(infile_path, min_seq_length):
    handle_elts = infile_path.split(".")
    handle_mod = "len" + str(args.min_seq_length)
    handle_elts.insert(-1,handle_mod)
    outfile_path = ".".join(handle_elts)
    return outfile_path

# parse incoming argument
parser = argparse.ArgumentParser()
parser.add_argument("fasta_file", help="Sequences in FASTA format")
parser.add_argument("min_seq_length", help="Minimum sequence length (as an integer) to keep", type=int)
args = parser.parse_args()

# open up the fasta
fasta_file = open(args.fasta_file, 'r')
output_fasta_handle = build_output_handle(args.fasta_file, args.min_seq_length)
print output_fasta_handle
output_fasta = open(output_fasta_handle, 'w')

seqs_of_min_len = 0
# extract seqs by length
for seq_record in SeqIO.parse(fasta_file, "fasta"):
    if len(seq_record.seq) >= args.min_seq_length:
        SeqIO.write(seq_record, output_fasta, "fasta")
        seqs_of_min_len += 1

print seqs_of_min_len, "sequences of at least length", args.min_seq_length, "were written to", output_fasta_handle
fasta_file.close()
output_fasta.close()
