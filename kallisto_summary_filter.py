#!/usr/bin/env python

# kallisto_summary_minmax.py

"""Print the max values and LQ value for est_counts and num_samples_found
in a kallisto summary file from count_kallisto_tpm_combined.py

Also summary metrics for the whole set:
# 3. Quartile statistics for norm counts & # samples

"""

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="summary csv_file", type=str)
parser.add_argument("-c", "--counts_file", help="norm-factored counts csv_file that matches summary file", type=str)
parser.add_argument("-l", "--min_count", help="Minimum number of samples this contig appears in, suggest at least 8", type=int)
parser.add_argument("-e", "--min_exp", help="Minimum average (normalized) expression", type=int)
parser.add_argument("-o", "--out_file", help="Output csv_file with rows from counts_file passing input criteria", type=str)
args = parser.parse_args()

# open files and collect variables:
input_csv = open(args.file, 'r')
counts_file = open(args.counts_file, 'r')
output_csv = open(args.out_file, 'w')
min_exp = int(args.min_exp)
min_count = int(args.min_count)

headers = input_csv.readline() # skip header
counts_header = counts_file.readline()
# "contig_id,avg_count,num_samples_found"
# turn those headers back around!
output_csv.write(counts_header)

for line in input_csv:
	counts_line = counts_file.readline()
	line_elts = line.split(",")
	contig_id = line_elts[0]
	avg_expression = float(line_elts[1])
	num_samples_found = int(line_elts[2].strip())
	if avg_expression >= min_exp and num_samples_found >= min_count:
		# counts_elts = counts_line.split(",")
		# count_id = counts_elts[0]
		count_id = counts_line.split(",")[0]
		assert count_id == contig_id
		output_csv.write(counts_line)
		# output_csv.write(line)
