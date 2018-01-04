#!/usr/bin/env python

# count_kallisto_tpm_combined.py

# we want to output two metrics:
# 1. mean tpm
# 2. no. of samples found in (or % of samples found in)

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="tsv_file", type=str)
parser.add_argument("-o", "--output_csv", help="output_csv", type=str)
args = parser.parse_args()

input_tsv = open(args.file, 'r')

output_csv = open(args.output_csv, 'w')
output_headers = "contig_id,mean_tpm,no_samples\n"
output_csv.write(output_headers)

headers = input_tsv.readline()
for line in input_tsv:
	# reset values
	total_tpm = 0.0
	mean_tpm = 0
	pos_samples = 0
	# collect basic info
	line_elts = line.split("\t")
	contig_id = line_elts[0]
	num_samples = (len(line_elts) - 1) # number of samples is number of cols - 1 (the ids)
	# go through each sample and count:
	for count in line_elts[1:]:
		count = float(count.strip())
		if count > 0:
			total_tpm += count
			pos_samples += 1
	mean_tpm = total_tpm / num_samples

	output_csv.write(contig_id + "," + str(mean_tpm) + "," + str(pos_samples) + "\n")
