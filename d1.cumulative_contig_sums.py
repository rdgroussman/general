#!/usr/bin/env python

# d1.cumulative_count_sums.py

"""
This script will calculate the cumulative, nested contig sums for the
combined counts file output by taxify_dmnd_lca_taxids.py

We will leverage nested taxonomic information to determine the cumulative sums
for each tax_id inclusive of its children tax_ids

Input file structure example: /Users/rgroussman/data/SCOPE/diel1/assemblies/diamond/d1.all_times.bf100.id99.vs_MarRef.uniq_taxids.taxified.csv

tax_id,tax_name,rank,count
267565,Skeletonema grethae,species,36,
267566,Skeletonema dohrnii,species,133,
267567,Skeletonema marinoi,species,376,
905079,Guillardia theta CCMP2712,below_species,1084,


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
ranks = ["root","below_root","superkingdom","below_superkingdom","below_below_superkingdom","below_below_below_superkingdom","below_below_below_below_superkingdom","kingdom","below_kingdom","below_below_kingdom","below_below_below_kingdom","below_below_below_below_kingdom","below_below_below_below_below_kingdom","subkingdom","phylum","below_phylum","below_below_phylum","below_below_below_phylum","below_below_below_below_phylum","below_below_below_below_below_phylum","below_below_below_below_below_below_phylum","below_below_below_below_below_below_below_phylum","below_below_below_below_below_below_below_below_phylum","below_below_below_below_below_below_below_below_below_phylum","below_below_below_below_below_below_below_below_below_below_phylum","subphylum","below_subphylum","below_below_subphylum","below_below_below_subphylum","below_below_below_below_subphylum","below_below_below_below_below_subphylum","below_below_below_below_below_below_subphylum","below_below_below_below_below_below_below_subphylum","below_below_below_below_below_below_below_below_subphylum","superclass","class","below_class","below_below_class","below_below_below_class","subclass","below_subclass","infraclass","below_infraclass","below_below_infraclass","below_below_below_infraclass","below_below_below_below_infraclass","superorder","below_superorder","order","below_order","below_below_order","suborder","superfamily","family","below_family","below_below_family","subfamily","tribe","genus","below_genus","species_group","species","below_species","below_below_species","subspecies","below_subspecies","varietas","below_varietas","forma"]
# now we'll go through each tax_id counts row:
for row in counts_data:
	tax_id = row["tax_id"]
	# for each 'sample_id', we will add the count to the tax_id in the field for any of the ranks:
	count = int(row["count"].strip())
	# and we'll add that count to each of the parent ranks (and the identity rank)
	for rank in ranks:
		if tax_id in TaxaDict:
			rank_id = TaxaDict[tax_id][rank]
		else:
			rank_id = ""
		# here is where we increment the count:
		if len(rank_id) > 0:
			CountsDict[rank_id]["total_counts"] += count

# output:
output_csv = open(args.output, 'w')
# this is the same ordering as sample_ids above, but with tax_id, and S6C1 -> S06C1 fixed, etc
header = ",".join(["tax_id","rank","tax_name","total_counts"])
output_csv.write(header + "\n")

for tax_id in CountsDict.keys():
	outline_elts = [str(tax_id)]
	outline_elts.append(str(CountsDict[tax_id]["rank"]))
	outline_elts.append(str(CountsDict[tax_id]["tax_name"]))
	outline_elts.append(str(CountsDict[tax_id]["total_counts"]))
	output_csv.write(",".join(outline_elts) + "\n")
