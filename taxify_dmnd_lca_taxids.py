#!/usr/bin/env python

# taxify_taxa_fracts.py

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-t", "--taxa_csv", help="Input taxa.csv", type=str)
parser.add_argument("-i", "--input", help="Input uniq_c.txt", type=str)
parser.add_argument("-o", "--out_file", help="Specify the name of the outfile", type=str)
args = parser.parse_args()

TAXA_CSV_PATH = args.taxa_csv
TAXA_COUNTS_PATH = args.input
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
# write the outfile header:
out_csv.write("tax_id,tax_name,rank,count\n")

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

# open the taxid counts file and iterate through each line:
taxid_counts = open(TAXA_COUNTS_PATH, 'r')
for line in taxid_counts:
	line_elts = line.split()
	count = line_elts[0]
	tax_id = line_elts[1].strip()
	if tax_id in TaxDict.keys():
		tax_name = TaxDict[tax_id]["tax_name"]
		rank = TaxDict[tax_id]["rank"]
	elif tax_id not in TaxDict.keys():
		tax_name = "Unk"
		rank = "Unk"
	# write out in this order: "tax_id,tax_name,rank,count"
	outline = ",".join([str(tax_id),tax_name,rank,str(count),"\n"])
	out_csv.write(outline)

# now open taxa_fracts:
# with open(TAXA_FRACTS_PATH, 'r') as taxa_fracts:
# 	for line in taxa_fracts:
# 		tax_id = line.strip()
# 		added_values = []
# 		# column_names =
# 		for column in ["tax_name"]:
# 			if tax_id in TaxDict:
# 				added_values.append(TaxDict[tax_id][column])
# 			else:
# 				TaxDict[tax_id] = {}
# 				TaxDict[tax_id]["tax_name"] = "NA"
# 				TaxDict[tax_id]["tax_id"] = tax_id
# 				TaxDict[tax_id]["rank"] = "NA"
# 				for i in range(len(column_data)):
# 					TaxDict[tax_id][column_data[i][1].strip('""')] = "NA"
# 				added_values.append(TaxDict[tax_id][column])
# 		for column in kept_fields:
# 			group_tax = TaxDict[tax_id][column]
# 			if group_tax != '':
# 				if group_tax in TaxDict:
# 					added_values.append(TaxDict[group_tax]["tax_name"])
# 				else:
# 					added_values.append("NA")
# 		merged_columns = ",".join(added_values)
# 		out_csv.write(tax_id + "," + merged_columns + "\n")
