#!/usr/bin/env python

# rank_counts_from_pplacer.py

""""Takes in guppy csv output from pplacer, and will return the counts for
each RANK. e.g., how many classifications are placed in phylum, class, order, etc """


import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="csv_file", type=str)
parser.add_argument("-t", "--taxa_csv", help="taxa_csv file", type=str)
parser.add_argument("-o", "--output_csv", help="output_summary_csv", type=str)
args = parser.parse_args()

guppy_csv = open(args.file, 'r')
taxa_csv = open(args.taxa_csv, 'r')
output_csv = open(args.output_csv, 'w')

# process taxa_csv into a dictionary:
TaxaDict = {}
rank_orders = taxa_csv.readline() # skip headers:
for line in taxa_csv:
	line_elts = line.split(",")
	tax_id = int(line_elts[0].strip('"'))
	rank = line_elts[2]
	TaxaDict[tax_id] = rank

# We'll store the counts for each rank in a dictionary
TaxaCounts = {}
# now go through the guppy csv file and increment counts for each rank:
"origin,name,multiplicity,edge_num,like_weight_ratio,post_prob,likelihood,marginal_like,distal_length,pendant_length,classification,map_ratio,map_overlap,map_identity"
headers = guppy_csv.readline()
for line in guppy_csv:
	line_elts = line.split(",")
	classification = int(line_elts[10].strip())
	if classification in TaxaDict.keys():
		rank = TaxaDict[classification]
		if rank in TaxaCounts.keys():
			TaxaCounts[rank] += 1
		elif rank not in TaxaCounts.keys():
			TaxaCounts[rank] = 1

# now we output the ranks in order given from taxa.csv:
for rank in rank_orders.split(",")[4:]:
	if rank in TaxaCounts.keys():
		outline = rank + "," + str(TaxaCounts[rank]) + "\n"
		output_csv.write(outline)
