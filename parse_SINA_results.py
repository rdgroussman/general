#!/usr/bin/env python


# parse_SINA_results.py

import argparse
import csv

parser = argparse.ArgumentParser()
parser.add_argument("input_csv", help="SINA-formatted csv delimited file")
parser.add_argument("-o", "--output_handle", help="Output handle", type=str)
parser.add_argument("-r", "--ranked_output", help="Output counts for all OTUs in each of the six ranks", action="store_true")
parser.add_argument("-k", "--krona_output", help="Output a format suitable for input into Krona", action="store_true")
args = parser.parse_args()


def parse_taxonomy(taxonomy):
	"""
	Parses semicolon-separated taxonomy values.
	The taxonomic paths are standardized to six ranks; Domain, Phylum, Class, Order, Family and Genus:

	Extremely limited example:
		Eukaryota;SAR;Alveolata;
		Eukaryota;
		Eukaryota;
		Eukaryota;SAR;Alveolata;Protalveolata;
		Eukaryota;SAR;Alveolata;Protalveolata;Syndiniales;
		Eukaryota;SAR;Alveolata;Protalveolata;Syndiniales;Syndiniales Group II;
	"""

	tax_elts = taxonomy.split(";")
	tax_len = len(tax_elts)
	for i in range(tax_len):
		if i == 0:
			add_to_dict('domain', tax_elts[i])
		if i == 1:
			add_to_dict('phylum', tax_elts[i])
		if i == 2:
			add_to_dict('class', tax_elts[i])
		if i == 3:
			add_to_dict('order', tax_elts[i])
		if i == 4:
			add_to_dict('family', tax_elts[i])
		if i == 5:
			add_to_dict('genus', tax_elts[i])

def add_to_dict(rank, name):
	"""Add the names found for the different ranks and increment a tally if they're found more than once"""

	if name not in TaxDict[rank].keys():
		TaxDict[rank][name] = 1
	elif name in TaxDict[rank].keys():
		TaxDict[rank][name] += 1

def catalog_taxonomy(taxonomy):
	"""
	Store the whole taxonomy and the counts for each:
	e.g.:
	'Eukaryota;SAR;Alveolata;Protalveolata;Syndiniales;Syndiniales Group II':5
	"""

	if taxonomy not in KronaTaxDict.keys():
		KronaTaxDict[taxonomy] = 1
	elif taxonomy in KronaTaxDict.keys():
		KronaTaxDict[taxonomy] += 1


# Initialize the taxonomy dictionary
TaxDict = {'domain':{}, 'phylum':{}, 'class':{}, 'order':{}, 'family':{}, 'genus':{}}
KronaTaxDict = {} # for use with args.krona_output == True
# "align_bp_score_slv,align_cutoff_head_slv,align_cutoff_tail_slv,align_family_slv,align_filter_slv,align_ident_slv,align_log_slv,align_quality_slv,align_startpos_slv,align_stoppos_slv,aligned_slv,full_name,lca_tax_slv,nearest_slv,nuc,nuc_gene_slv,turn_slv"
input_csv = open(args.input_csv, 'r')
header = input_csv.readline() # skip the header
# this is necessary to avoid the commas within double-quotes
for line in csv.reader(input_csv, quotechar='"', delimiter=",", quoting=csv.QUOTE_ALL, skipinitialspace=True):
	taxonomy = line[12].strip(";") # nested taxonomy is given in field 12 with semicolon separator, e.g: "Eukaryota;SAR;Alveolata;Protalveolata;Syndiniales;Syndiniales Group II;"
	parse_taxonomy(taxonomy) # sends this to the taxonomy parser
	if args.krona_output == True:
		catalog_taxonomy(taxonomy)

# Output the results!
if args.ranked_output == True:
	rank_ordering = ['domain', 'phylum', 'class', 'order', 'family', 'genus']
	for rank in rank_ordering:
		outfile_path = args.output_handle + "." + rank + ".csv"
		outfile = open(outfile_path, 'w')
		for taxon in TaxDict[rank].keys():
			outline = taxon + "," + str(TaxDict[rank][taxon]) + "\n"
			outfile.write(outline)
		outfile.close()

if args.krona_output == True:
	# output in krona-friendly format:
	outfile_path = args.output_handle + ".krona.tsv"
	outfile = open(outfile_path, 'w')
	for taxonomy in sorted(KronaTaxDict.keys()):
		outline = str(KronaTaxDict[taxonomy]) + "\t" + taxonomy.replace(";", "\t") + "\n"
		outfile.write(outline)
	outfile.close()
