#!/usr/bin/env python

# taxify_taxa_fracts.py

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-t", "--taxa_csv", help="Input taxa.csv", type=str)
parser.add_argument("-i", "--input", help="Input taxa_fracts.csv", type=str)
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

kept_fields = ["phylum","class"]


def parse_taxacsv_header(header):

	column_data = []
	header_elts = header.split(",")
	for i in range(len(header_elts)):
		elt = header_elts[i].strip('"')
		if elt in kept_fields:
			column_data.append((i, elt))
	return column_data

TaxDict = {}
with open(args.taxa_csv, 'r') as taxa_csv:
	header = next(taxa_csv)
	column_data = parse_taxacsv_header(header)
	for line in taxa_csv:
		line_elts = line.split(",")
		tax_id = line_elts[0].strip('""')
		rank = line_elts[2].strip('""')
		tax_name = line_elts[3].strip('""')
		TaxDict[tax_id] = {}
		TaxDict[tax_id]["tax_name"] = tax_name
		TaxDict[tax_id]["tax_id"] = tax_id
		TaxDict[tax_id]["rank"] = rank
		for i in range(len(column_data)):
			TaxDict[tax_id][column_data[i][1].strip('""')] = line_elts[int(column_data[i][0])].strip('""')

# now open taxa_fracts:
with open(TAXA_FRACTS_PATH, 'r') as taxa_fracts:
	for line in taxa_fracts:
		tax_id = line.strip()
		added_values = []
		# column_names =
		for column in ["tax_name"]:
			added_values.append(TaxDict[tax_id][column])
		for column in kept_fields:
			group_tax = TaxDict[tax_id][column]
			if group_tax != '':
				added_values.append(TaxDict[group_tax]["tax_name"])
		merged_columns = ",".join(added_values)
		out_csv.write(tax_id + "," + merged_columns + "\n")
