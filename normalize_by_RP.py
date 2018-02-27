#!/usr/bin/env python

# normalize_by_RP.py

"""
This script is intended to normalize the per-group counts produced
by a script like count_pplacer_csv_by_taxonomy.py by similar
counts from ribosomal proteins or housekeeping genes.

The input csvs should have identical sample names but this script
will normalize groups (e.g, Bacillariphyta, Archaea, etc) only when present in both.
E.g, the counts to be normalized can have a column for 'Virus' counts while the
normalizing factor does not need to.

Outputs a CSV file with normalized counts for groups present in both counts files.

"""

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-n", "--numerator", help="Path to numerator csv", type=str)
parser.add_argument("-d", "--denominator", help="Path to denominator csv", type=str)
parser.add_argument("-o", "--output", help="Name of output csv", type=str)

# Parse the arguments given at the command line
args = parser.parse_args()

# load the numerator counts file
numerator_counts = open(args.numerator, 'r')

# capture the header:
header = numerator_counts.readline().strip()
# get the set of groups in numerator counts
num_header_elts = header.split(",")
tax_groups = num_header_elts[1:]

# build a nested dictionary of counts, with each sample as its own dictionary:
NumeratorDict = {}
for line in numerator_counts:
	line_elts = line.strip().split(",")
	sample = line_elts[0]
	NumeratorDict[sample] = {}
	# for each sample dictionary, give the counts for the associated group:
	for i in range(len(tax_groups)):
		NumeratorDict[sample][tax_groups[i]] = float(line_elts[i+1])

# load the denominator counts file
denominator_counts = open(args.denominator, 'r')
header = denominator_counts.readline().strip()
# get the set of groups in numerator counts
den_header_elts = header.split(",")
den_tax_groups = den_header_elts[1:]

# build a nested dictionary of counts, with each sample as its own dictionary:
DenominatorDict = {}
for line in denominator_counts:
	line_elts = line.strip().split(",")
	sample = line_elts[0]
	DenominatorDict[sample] = {}
	# for each sample dictionary, give the counts for the associated group:
	for i in range(len(tax_groups)):
		DenominatorDict[sample][tax_groups[i]] = float(line_elts[i+1])

# look at the intersection to determine the groups to perform normalization on:
# shared taxonomic groups:
common_tax = set(tax_groups) & set(den_tax_groups)
# shared samples:
common_samples = set(NumeratorDict.keys()) & set(DenominatorDict.keys())

# now we go through and divide the values in the NumeratorDict by the DenominatorDict
# storing the new values in NormDict and giving the value 'NaN' when the value in DenominatorDict is 0
NormDict = {}
for sample in common_samples:
	NormDict[sample] = {}
	for tax in common_tax:
		if NumeratorDict[sample][tax] > 0.0 and DenominatorDict[sample][tax] == 0.0:
			NormDict[sample][tax] = "NaN"
		elif NumeratorDict[sample][tax] == 0.0 and DenominatorDict[sample][tax] == 0.0:
			NormDict[sample][tax] = 0
		else:
			NormDict[sample][tax] = NumeratorDict[sample][tax] / DenominatorDict[sample][tax]


# finally, output the results:
outfile = open(args.output, 'w')
outheader = "sample_name," + ",".join(sorted(common_tax)) + ",\n"
outfile.write(outheader)
for sample in sorted(common_samples):
	outlist = [sample]
	for tax in sorted(common_tax):
		outlist.append(str(NormDict[sample][tax]))
	outline = ",".join(outlist) + ",\n"
	outfile.write(outline)
