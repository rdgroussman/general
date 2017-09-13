#!/usr/bin/env python


# parse_polgyenic_presence.py

"""
# GIVEN A LIST OF GENE IDS
For each gene in the list, the script will comb through a series of FASTA files
and maintain a set of tax_ids that are present in each.

list: txt file with sequential names of FASTA files:
	uS3/uS3.MarRef2Plus.query.fasta
	uS8/uS8.MarRef2Plus.query.fasta
	uS10/uS10.MarRef2Plus.query.fasta
	uS17/uS17.MarRef2Plus.query.fasta
	uS19/uS19.MarRef2Plus.query.fasta
	uL2/uL2.MarRef2Plus.query.fasta

* cool trick: This script will work on the deflines only, doesn't need whole fasta files.

"""


import argparse
import re

# from Bio import Phylo

parser = argparse.ArgumentParser()
parser.add_argument("input_list", help="Txt file list of input FASTAs")
parser.add_argument("-o", "--out_file", help="Output file name", type=str)
args = parser.parse_args()

# here's our regular expression for finding taxids:
tax_pattern = re.compile("(tax[0-9]+)")

def process_fasta(fasta):

	# Open the file
	fasta_file = open(fasta, 'r')

	# Iterate through lines, only read deflines (starts with >)
	for line in fasta_file:
		if line.startswith(">"):
		    line_elts = line.split("_")
		    for elt in line_elts:
		        if tax_pattern.match(elt):
		            tax_id = tax_pattern.findall(elt)[0][3:]
		            taxid_dict[fasta].add(tax_id)

	# close the fasta!
	fasta_file.close()

# Open up and parse the list
listfile = open(args.input_list, 'r')
outfile = open(args.out_file, 'w')
fasta_list = []
taxid_dict = {}
for line in listfile:
	fasta = line.strip()
	fasta_list.append(fasta)
	taxid_dict[fasta] = set([])

# Go through each FASTA file in the list:
for fasta in fasta_list:
	process_fasta(fasta)

unique_set = taxid_dict[fasta_list[0]].copy() # unique set starts as a copy of the first set (otherwise no intersection!)
global_set = set([]) # empty set
# Report results:
for fasta in fasta_list:
	# Print the len(set) for each file
	print fasta + "contains " + str(len(taxid_dict[fasta])) + " taxa"
	# add it to the intersection set:
	unique_set = unique_set.intersection(taxid_dict[fasta])
	global_set = global_set.union(taxid_dict[fasta])

print "Global set contains a total of " + str(len(global_set)) + " taxa"
print "Unique set contains a total of " + str(len(unique_set)) + " taxa"

print "Unique taxids written out to " + args.out_file
for i in unique_set:
	outfile.write(str(i) + "\n")
outfile.close()
