#!/usr/bin/env python


# taxify_taxa_fracts.py

TAXA_CSV_PATH="/Users/rgroussman/data/SCOPE/diel1/single_gene/DGAT1/taxa.csv"
TAXA_FRACTS_PATH="/Users/rgroussman/data/SCOPE/diel1/single_gene/DGAT1/summed_csv/taxa_fracts.csv"
OUT_CSV_PATH="/Users/rgroussman/data/SCOPE/diel1/single_gene/DGAT1/summed_csv/taxified_taxa_fracts.csv"

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
		TaxDict[tax_id] = [ tax_name, rank ]
tax_id = None
#
# print TaxDict
# tax_id = '31353'
# print TaxDict[tax_id]

# now open taxa_fracts:
with open(TAXA_FRACTS_PATH, 'r') as taxa_fracts:
	next(taxa_fracts)
	for line in taxa_fracts:
		tax_id = str(line.split(',')[2])
		added_columns = ",".join(TaxDict[tax_id])
		out_csv.write(line.strip() + "," + added_columns + "\n")
