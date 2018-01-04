#!/usr/bin/env python


# kallisto_norm_factors.py

# takes in a big file of counts as well as a table with norm factors
# outputs a table of counts modified by norm factors

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-n", "--normfactors_csv", help="Specify path to an csv file with norm factors for your samples.", type=str)
parser.add_argument("-o", "--outfile", help="Specify name of output file.", type=str)
parser.add_argument("input_tab", help="tab-formatted input file")
args = parser.parse_args()

def fix_sample_name(sample_name):
	"""Fixes some sample names to standardize format.
	E.g. S8C1 > S08C1, etc """

	if sample_name.startswith("S6C1"):
		sample_name = sample_name.replace("S6C1", "S06C1")
	if sample_name.startswith("S7C1"):
		sample_name = sample_name.replace("S7C1", "S07C1")
	if sample_name.startswith("S8C1"):
		sample_name = sample_name.replace("S8C1", "S08C1")
	return sample_name

def initialize_norm_counts_dict():

	NormFactorsDict = {}
	normfactors_csv = open(args.normfactors_csv, 'r')
	for line in normfactors_csv:
		line_elts = line.split(",")
		sample = line_elts[0]
		norm_factor = float(line_elts[1].strip())
		NormFactorsDict[sample] = norm_factor
	return NormFactorsDict

def get_normfactor(i):
	sample = SampleOrder[i]
	norm_factor = NormFactorsDict[sample]
	return norm_factor

# Parse the normfactors_csv:
NormFactorsDict = initialize_norm_counts_dict()

# Now go through the target tab file line by line and modify the counts:
input_counts = open(args.input_tab, 'r')
outfile = open(args.outfile, 'w')

# we'll use the header to determine the right ordering of the samples. Make a dictionary:
SampleOrder = {}
header = input_counts.readline()
header_elts = header.split("\t")
for i in range(len(header_elts)):
	SampleOrder[i] = header_elts[i].strip()
outfile.write(",".join(header_elts))

# now we'll iterate through the rest of the lines in input_counts
for line in input_counts:
	line_elts = line.split("\t")
	contig_id = line_elts[0]
	out_line = [contig_id]
	for i in range(len(line_elts)):
		if i > 0:
			norm_factor = get_normfactor(i)
			corrected_count = norm_factor * float(line_elts[i].strip())
			out_line.append(str(corrected_count))
	out_line.append("\n")
	outfile.write(",".join(out_line))
