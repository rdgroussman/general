#!/usr/bin/env python

# check_hmmtab_for_dupes.py
"""Goes through the input hmm_out.tab file, and checks for the presence of
duplicates. If it finds a duplicate, only write out the first incident.
"""

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("input_tab", help="Specify input table from hmmsearch")
args = parser.parse_args()


unique_hits = {}
check_hmm_table = open(args.input_tab, 'r')
for line in check_hmm_table:
	if line[0] != '#':
		hit_id = line.split('-')[0].strip()
		if hit_id not in unique_hits:
			unique_hits[hit_id] = line
		elif hit_id in unique_hits:
			print "DUPLICATE HIT FOUND IN ", args.input_tab
			print "1:", unique_hits[hit_id]
			print "2:", line

check_hmm_table.close()
# /mnt/nfs/ryan/scripts/check_hmmtab_for_dupes.py
