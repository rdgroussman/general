#!/usr/bin/env python


# filter_guppy_csv_by_pendant_length.py

# A simple script. Its purpose? Take in a 'guppy' csv and a pendant length cutoff.
# Output the lines with a pendant length below a certain cutoff threshold.

# Pendant length is generally > 0.7 (to a maximum of 2)

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="Input guppy csv from pplacer", type=str)
parser.add_argument("-o", "--out_file", help="Specify the name of the outfile (csv)", type=str)
parser.add_argument("-c", "--cutoff", help="Specify the pendant length cutoff", type=float)
args = parser.parse_args()

out_csv = open(args.out_file,'w')
cutoff = args.cutoff
line_counter = 0
discarded_counter = 0

taxa_csv = open(args.input, 'r')
for line in taxa_csv:
	line_counter += 1 # increment line counter
	if line_counter == 1:
		out_csv.write(line) # write out the header
	elif line_counter > 1:
		line_elts = line.split(",")
		pendant_length = float(line_elts[9])
		if pendant_length <= cutoff:
			out_csv.write(line)
		elif pendant_length > cutoff:
			discarded_counter += 1

print "Filtering", args.input, "by pendant length", args.cutoff
print "Removed", discarded_counter, "out of", (line_counter -1), "reads"
print "Saved filtered output to", args.out_file
out_csv.close()
