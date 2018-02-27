#!/usr/bin/env python

# count_kallisto_tpm_combined.py

# we want to output two metrics:
# 1. mean tpm
# 2. no. of samples found in (or % of samples found in)

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="tsv_file", type=str)
parser.add_argument("-o", "--output_csv", help="output_summary_csv", type=str)
parser.add_argument("-t", "--output_tsv", help="output_filtered_tsv", type=str)
parser.add_argument("-l", "--min_count", help="Minimum number of samples this contig appears in, suggest at least 8", type=int)
args = parser.parse_args()

"""
Script should do these thing:
# filter out all contigs with 32 or more zeros (or given by --min_count)
# record how many contigs have 0 hits
# and how many contigs are cut with the 32 threshold (or given by --min_count)
"""

# Initialize and load input and output files:
input_tsv = open(args.file, 'r')

output_tsv = open(args.output_tsv, 'w')

output_csv = open(args.output_csv, 'w')
output_headers = "contig_id,est_counts,num_samples_found\n"
output_csv.write(output_headers)


# counters for all-zeroes and 32+ (or whatever interger given) zeroes:
allzero_counter = 0
below_min_counter = 0
total_contigs = 0
min_count = args.min_count

headers = input_tsv.readline()
output_tsv.write(headers)
for line in input_tsv:
	total_contigs += 1
	# reset values
	total_counts = 0.0
	mean_counts = 0
	pos_samples = 0
	# collect basic info
	line_elts = line.split(",")
	contig_id = line_elts[0]
	num_samples = (len(line_elts) - 1) # number of samples is number of cols - 1 (the ids)
	# go through each sample and count:
	for count in line_elts[1:]:
		count = float(count.strip())
		if count > 0:
			total_counts += count
			pos_samples += 1
	mean_counts = total_counts / num_samples

	# increment counters and write out the line if it passes minimum presence criteria:
	if pos_samples == 0:
		allzero_counter += 1
	elif pos_samples < min_count:
		below_min_counter += 1
	elif pos_samples >= min_count:
		output_tsv.write(line)

	# write out to the summary_csv:
	output_csv.write(contig_id + "," + str(mean_counts) + "," + str(pos_samples) + "\n")

# print some summary statistics out to screen:
print "Total of ", str(total_contigs), "in", args.file
print "Number of contigs with zero-counts in all samples: ", str(allzero_counter)
print "Number of contigs below the minimum presence threshold of", str(min_count), ":", str(below_min_counter)


# close files:
input_tsv.close()
output_tsv.close()
output_csv.close()
