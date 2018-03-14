#!/usr/bin/env python

# count_kallisto_tpm_combined.py

# we want to output two metrics for each contig:
# 1. mean tpm or est_counts (given upstream by a cat $dir/abundance.tsv | awk -F' ' command)
# 2. no. of samples found in (or % of samples found in)

import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="tsv_file", type=str)
parser.add_argument("-o", "--output_csv", help="output_summary_csv", type=str)
parser.add_argument("-t", "--summary_text", help="summary_text", type=str)
parser.add_argument("-l", "--min_count", help="Minimum number of samples this contig appears in, suggest at least 8", type=int)
parser.add_argument("-e", "--min_exp", help="Minimum average expression for this contig across all samples (e.g., lowest quartile value)", type=int)
args = parser.parse_args()

"""
Script should do these thing:
# we want to output two metrics for each contig:
# 1. mean tpm or est_counts (given upstream by a cat $dir/abundance.tsv | awk -F' ' command)
# 2. no. of samples found in (or pct of samples found in)

And one summary metric:
# 3. Quartile statistics for norm counts & # samples
"""

# Initialize and load input and output files:
input_tsv = open(args.file, 'r')

output_txt = open(args.summary_text, 'w')

output_csv = open(args.output_csv, 'w')
output_headers = "contig_id,avg_count,num_samples_found\n"
output_csv.write(output_headers)

# we will store our collected average counts in this list:
avg_count_list = []

# counters for all-zeroes and 32+ (or whatever interger given) zeroes:
allzero_counter = 0
below_min_counter = 0
total_contigs = 0
min_count = args.min_count

headers = input_tsv.readline() # skip the header
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
	avg_count_list.append(mean_counts)

	# increment counters and write out the line if it passes minimum presence criteria:
	if pos_samples == 0:
		allzero_counter += 1
	elif pos_samples < min_count:
		below_min_counter += 1
	# elif pos_samples >= min_count:
		# output_tsv.write(line)

	# write out to the summary_csv:
	output_csv.write(contig_id + "," + str(int(round(mean_counts))) + "," + str(pos_samples) + "\n")

# print some summary statistics to the summary text:
outstring = "Total of " + str(total_contigs) + " in " + args.file + "\n"
output_txt.write(outstring)
outstring = "Number of contigs with zero-counts in all samples: " + str(allzero_counter) + "\n"
output_txt.write(outstring)
outstring = "Number of contigs below the minimum presence threshold of " + str(min_count) + ": " + str(below_min_counter) + "\n"
output_txt.write(outstring)
outstring = "####### Percentile statistics #######" + "\n"
output_txt.write(outstring)
outstring = "Minimum average expression: " + str(np.percentile(avg_count_list,0)) + "\n"
output_txt.write(outstring)
outstring = "1st quartile average expression: " + str(np.percentile(avg_count_list,25)) + "\n"
output_txt.write(outstring)
outstring = "2nd quartile (median) average expression: " + str(np.percentile(avg_count_list,50)) + "\n"
output_txt.write(outstring)
outstring = "3rd quartile average expression: " + str(np.percentile(avg_count_list,75)) + "\n"
output_txt.write(outstring)
outstring = "Maximum average expression: " + str(np.percentile(avg_count_list,100)) + "\n"
output_txt.write(outstring)

# close files:
input_tsv.close()
output_txt.close()
output_csv.close()
