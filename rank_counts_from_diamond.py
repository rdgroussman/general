#!/usr/bin/env python

# rank_counts_from_diamond.py

""""Takes in LCA results from DIAMOND, and will return the counts for
each RANK. e.g., how many classifications are placed in phylum, class, order, etc """


import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="tab_file", type=str)
parser.add_argument("-t", "--taxa_csv", help="taxa_csv file", type=str)
parser.add_argument("-o", "--output_csv", help="output_summary_csv", type=str)
args = parser.parse_args()

diamond_lca = open(args.file, 'r')
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
# and we'll add a special taxid for 0, the non-placement result:
TaxaDict[int("0")] = "Unclassified"

# We'll store the counts for each rank in a dictionary
TaxaCounts = {}
# now go through the guppy csv file and increment counts for each rank:
# "contig_id,tax_id,evalue"
for line in diamond_lca:
	line_elts = line.split()
	classification = int(line_elts[1].strip())
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
# and a special case for the unclassified:
if "Unclassified" in TaxaCounts.keys():
	outline = "Unclassified" + "," + str(TaxaCounts['Unclassified']) + "\n"
	output_csv.write(outline)
