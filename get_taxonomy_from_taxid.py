#!/usr/bin/env python

# taxify_taxa_fracts.py

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-t", "--taxa_csv", help="Input taxa.csv", type=str)
parser.add_argument("-i", "--input", help="Input tax_id", type=str)
parser.add_argument("-o", "--out_file", help="Specify the name of the outfile", type=str)
args = parser.parse_args()

TAXA_CSV_PATH = args.taxa_csv
TAXA_FRACTS_PATH = args.input
OUT_CSV_PATH = args.out_file

# old hard-coded values
# TAXA_CSV_PATH="/Users/rgroussman/data/SCOPE/diel1/single_gene/DGAT1/taxa.csv"
# TAXA_FRACTS_PATH="/Users/rgroussman/data/SCOPE/diel1/single_gene/DGAT1/summed_csv/taxa_fracts.csv"
# OUT_CSV_PATH="/Users/rgroussman/data/SCOPE/diel1/single_gene/DGAT1/summed_csv/taxified_taxa_fracts.csv"

# get tax_id, rank, and tax_name from taxa.csv
# put latter two in dictionary under key tax_id

# declare the dictionary:
TaxDict = {}
out_csv = open(OUT_CSV_PATH,'w')


with open(TAXA_CSV_PATH, 'r') as taxa_csv:
	next(taxa_csv)
	for line in taxa_csv:
		line_elts = line.split(",")
		tax_id = line_elts[0].strip('""')
		rank = line_elts[2].strip('""')
		# get rid of the quotes around values in taxa.csv
		tax_name = line_elts[3].strip('""')
		phylum_id = line_elts[17].strip('""')
		class_id = line_elts[38].strip('""')
		TaxDict[tax_id] = [ tax_name, rank, phylum_id, class_id ]
tax_id = None

# print TaxDict
# tax_id = '31353'
# print TaxDict[tax_id]

# now open taxa_fracts:
with open(TAXA_FRACTS_PATH, 'r') as taxa_fracts:
	for line in taxa_fracts:
		tax_id = line.strip()
		phylum = ""
		if TaxDict[tax_id][2] != "":
			phylum_id = TaxDict[tax_id][2]
			phylum = TaxDict[phylum_id][0]
		class_ = ""
		if TaxDict[tax_id][3] != "":
			class_id = TaxDict[tax_id][3]
			class_ = TaxDict[class_id][0]

		added_columns = ",".join([phylum, class_])
		out_csv.write(line.strip() + "," + added_columns + "\n")
