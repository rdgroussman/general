#!/usr/bin/env python


# d1.cumulative_pfam_sums.py


"""
This script will calculate the cumulative, nested pfam sums for the
combined counts file output by count_pfam_by_taxid.py

We will leverage nested taxonomic information to determine the cumulative sums
for each tax_id inclusive of its children tax_ids

Input file structure example: /Users/rgroussman/data/SCOPE/diel1/assemblies/sql_db/data/d1.pfam_by_taxid.csv

tax_id,PF11081.7,PF14561.5,PF05161.12,... (7994 PFAMs)
267565,0,1,5,0,0, (etc)

This script will also go a little further and trim the pfam_id: PF11081.7 > PF11081, etc


"""

import argparse
import csv

parser = argparse.ArgumentParser()
parser.add_argument("-l", "--csv", help="List of *est_counts_by_tax_id.csv files", type=str)
parser.add_argument("-t", "--taxa_csv", help="NCBI taxa.csv file", type=str)
parser.add_argument("-o", "--output", help="Output CSV path", type=str)
args = parser.parse_args()

# go through and parse the taxa_csv file
TaxaDict = {}
CountsDict = {}
pfam_set = set([])

taxa_csv = csv.DictReader(open((args.taxa_csv),'r'))
# each tax_id gets its own dictionary from the row
for row in taxa_csv:
	tax_id = row["tax_id"]
	TaxaDict[tax_id] = row
	# also, populate the CountsDict:
	CountsDict[tax_id] = {"total_counts":0}
	CountsDict[tax_id]["rank"] = row["rank"]
	CountsDict[tax_id]["tax_name"] = row["tax_name"]

# then, load in the counts file and iterate through the taxids
counts_data = csv.DictReader(open((args.csv),'r'))
# sample_ids = ["total_counts","S11C1_A_1800","S11C1_C_1800","S14C1_B_2200","S14C1_C_2200","S15C1_B_200","S15C1_C_200","S16C1_A_600","S16C1_B_600","S17C1_A_1000","S17C1_B_1000","S18C1_A_1400","S18C1_C_1400","S19C1_A_1800","S19C1_C_1800","S20C1_B_2200","S20C1_C_2200","S21C1_B_200","S21C1_C_200","S22C1_A_600","S22C1_C_600","S23C1_B_1000","S23C1_C_1000","S24C1_A_1400","S24C1_B_1400","S26C1_A_1800","S26C1_C_1800","S28C1_B_2200","S28C1_C_2200","S29C1_A_200","S29C1_C_200","S30C1_A_600","S30C1_C_600","S31C1_A_1000","S31C1_C_1000","S32C1_B_1400","S32C1_C_1400","S33C1_A_1800","S33C1_C_1800","S34C1_B_2200","S34C1_C_2200","S35C1_A_200","S35C1_C_200","S06C1_A_600","S06C1_C_600","S07C1_A_1000","S07C1_B_1000","S08C1_B_1400","S08C1_C_1400"]

ranks = ["root","below_root","superkingdom","below_superkingdom","below_below_superkingdom","below_below_below_superkingdom","below_below_below_below_superkingdom","kingdom","below_kingdom","below_below_kingdom","below_below_below_kingdom","below_below_below_below_kingdom","below_below_below_below_below_kingdom","subkingdom","phylum","below_phylum","below_below_phylum","below_below_below_phylum","below_below_below_below_phylum","below_below_below_below_below_phylum","below_below_below_below_below_below_phylum","below_below_below_below_below_below_below_phylum","below_below_below_below_below_below_below_below_phylum","below_below_below_below_below_below_below_below_below_phylum","below_below_below_below_below_below_below_below_below_below_phylum","subphylum","below_subphylum","below_below_subphylum","below_below_below_subphylum","below_below_below_below_subphylum","below_below_below_below_below_subphylum","below_below_below_below_below_below_subphylum","below_below_below_below_below_below_below_subphylum","below_below_below_below_below_below_below_below_subphylum","superclass","class","below_class","below_below_class","below_below_below_class","subclass","below_subclass","infraclass","below_infraclass","below_below_infraclass","below_below_below_infraclass","below_below_below_below_infraclass","superorder","below_superorder","order","below_order","below_below_order","suborder","superfamily","family","below_family","below_below_family","subfamily","tribe","genus","below_genus","species_group","species","below_species","below_below_species","subspecies","below_subspecies","varietas","below_varietas","forma"]
# now we'll go through each tax_id counts row:
line_counter = 0
say_when = [10,25,50,100,200,300,400,500,600,700]

for row in counts_data:
	tax_id = row["tax_id"]
	line_counter += 1
	if line_counter in say_when:
		print line_counter
	# for each 'pfam_id', we will add the count to the tax_id in the field for any of the ranks:
	for pfam_id in row.keys():
		if pfam_id == "tax_id":
			pass
		else:
			# first, get the count:
			count = int(row[pfam_id].strip())
			# then trim the id:
			pfam_id = pfam_id.split(".")[0]
			# and we'll add that count to each of the parent ranks (and the identity rank)
			for rank in ranks:
				if tax_id in TaxaDict:
					rank_id = TaxaDict[tax_id][rank]
				else:
					rank_id = ""
				# here is where we increment the count:
				if len(rank_id) > 0:
					if pfam_id not in CountsDict[rank_id].keys():
						CountsDict[rank_id][pfam_id] = count
					else:
						CountsDict[rank_id][pfam_id] += count
			pfam_set.add(pfam_id)

# now, we'll make the pfam_set into a list so we can maintain ordering:
pfam_list = []
for pfam_id in pfam_set:
	pfam_list.append(pfam_id)

# now we'll open up our output file and build the header:
output_csv = open(args.output, 'w')
header = ["tax_id","rank","tax_name"]
for pfam_id in pfam_list:
	header.append(pfam_id)
header = ",".join(header) + "\n"
output_csv.write(header)

# now, go through each tax_id:
for tax_id in CountsDict.keys():
	outline = [str(tax_id)]
	outline.append(str(CountsDict[tax_id]["rank"]))
	outline.append(str(CountsDict[tax_id]["tax_name"]))
	for pfam_id in pfam_list:
		outline.append(str(CountsDict[tax_id].get(pfam_id,0)))
	outline = ",".join(outline) + "\n"
	output_csv.write(outline)
