#!/usr/bin/env python


# force_unique_headers.py

"""
This script takes a FASTA file as input.
Also takes additional arguments to use as sequence 'prefix'
In addition to adding the universal 'prefix', this script will make sure
that each header in the file is unique by referencing a unique set.
If a header is found to be a duplicate, it will append with _d1, _d2 etc
(for duplicate 1, duplicate 2, etc)

Practice example:

'input.fasta' has 341 sequences with the same header: ">HF7JC_8768_4_1"
A prefix is supplied with -p S28C1_B_2200_

Output fasta will contain the following sequence:
>S28C1_B_2200_HF7JC_8768_4_1
>S28C1_B_2200_HF7JC_8768_4_1_d1
>S28C1_B_2200_HF7JC_8768_4_1_d2
...
>S28C1_B_2200_HF7JC_8768_4_1_d340

"""


import argparse
parser = argparse.ArgumentParser()
parser.add_argument("input_fasta", help="FASTA-formatted input")
parser.add_argument("-p", "--prefix", help="Specify prefix for all sequences", type=str, default="")
parser.add_argument("-o", "--out_file", help="Output file name", type=str)
args = parser.parse_args()

def process_defline(line):

	# remove the pointy bracket '>'
	line = line[1:]
	# Check to see if the defline is in the set, make unique or add to set
	if line in HeaderSet:
		line = make_unique(line) # if it is, make it unique:
	elif line not in HeaderSet:
		HeaderSet.add(line)
	# Add bracket, prefix, and newline and return
	new_defline = ">" + defline_prefix + line + "\n"
	return new_defline

def make_unique(line):

	# UniqueDict keeps track of how many duplicates we've seen so far
	if line not in UniqueDict:
		UniqueDict[line] = 1
	elif line in UniqueDict:
		UniqueDict[line] += 1

	uniq_line = line + "_d" + str(UniqueDict[line])
	return uniq_line

in_file = open(args.input_fasta, 'r')
out_file = open(args.out_file, 'w')
defline_prefix = args.prefix

HeaderSet = set([])
UniqueDict = {}


for line in in_file:
	if line.startswith(">"):
		new_defline = process_defline(line.strip())
		out_file.write(new_defline)
	elif line.startswith(">") == False:
		out_file.write(line)
