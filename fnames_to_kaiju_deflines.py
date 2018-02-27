#!/usr/bin/env python

# fnames_to_kaiju_deflines.py

# Ryan Groussman
# Armbrust Lab, University of Washington
# January 2018

"""
This script converts the deflines from 'friendly names' format to
a format specified for input into kaiju custom database creation.

Input:

>SAR116_cluster_alpha_proteobacterium_HIMB100_tax909943_locID_2503127797_SAR116_HIMB100_00023710_seqID6
MKWGPGTLLFIHIAGSASEQMTEIPVAELVEGKGIKGDRYFNGTGNYSHIPDIREVTLIE
>SAR116_cluster_alpha_proteobacterium_HIMB100_tax909943_locID_2503127796_SAR116_HIMB100_00023700_seqID7
MLVSSPEDIATYICQIPKGGVVTPKKMRLDLARAKGADNSCPVSTGIFLRIAIEDVLRLF

Output style: counter_NCBItaxid
>1_1358
MAQQRRGGFKRRKKVDFIAANKIEVVDYKDTELLKRFISERGKILPRRVTGTSAKNQRKVVNAIKRARVMALLPFVAEDQN
>2_44689
MASTQNIVEEVQKMLDTYDTNKDGEITKAEAVEYFKGKKAFNPERSAIYLFQVYDKDNDGKITIKELAGDIDFDKALKEYKEKQAKSKQQEAEVEEDIEAFILRHNKDDNTDITKDELIQGFKETGAKDPEKSANFILTEMDTNKDGTITVKELRVYYQKVQKLLNPDQ
>3_352472
MKTKSSNNIKKIYYISSILVGIYLCWQIIIQIIFLMDNSIAILEAIGMVVFIS

"""

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("input_fasta", help="FASTA-formatted input file")
parser.add_argument("-o", "--output", help="Name of output file", type=str)
args = parser.parse_args()

input_fasta = open(args.input_fasta, 'r')
output_fasta = open(args.output, 'w')

def parse_header(defline):
	"""pick out only the NCBI taxid and return;
	given in the form '_tax#####'"""

	def_elts = defline.split("_")
	for elt in def_elts:
		if elt.startswith("tax"):
			taxid=elt.strip()[3:]
			if "/" in taxid: # some of the entries have a taxid followed by slash
				taxid = taxid.split("/")[0]
			return taxid

counter = 0
for line in input_fasta:
	if line.startswith(">"):
		counter += 1
		taxid = parse_header(line)
		newline = ">" + str(counter) + "_" + taxid + "\n"
		output_fasta.write(newline)
	else:
		output_fasta.write(line)
