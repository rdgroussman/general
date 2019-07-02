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
parser.add_argument("-l", "--csv", help="List of csv files", type=str)
parser.add_argument("-t", "--taxa_csv", help="NCBI taxa.csv file", type=str)
parser.add_argument("-o", "--output", help="Output CSV path", type=str)
args = parser.parse_args()


# we'll declare our discrete and cumulative groups to count here:

# discrete: we'll only count placements at this exact taxid:
"""
2759	Eukaryota	superkingdom
131567	cellular organisms	below_root
"""
discrete_groups = ["2759","131567"]


# we'll count placements in all of the nested taxa here:
"""
2864	Dinophyceae	class
33630	Alveolata	below_superkingdom
2830	Haptophyceae	below_superkingdom
33634	Stramenopiles	below_superkingdom
2836	Bacillariophyta	phylum
39119	Dictyochophyceae	class
33154	Opisthokonta	below_superkingdom
35675	Pelagophyceae	class
3041	Chlorophyta	phylum
33208	Metazoa	kingdom
543769	Rhizaria	below_superkingdom
431838	Intramacronucleata	subphylum
91989	Bolidophyceae	phylum
3027	Cryptophyta	class
2	Bacteria	superkingdom
5878	Ciliophora	below_below_superkingdom
10239	Viruses	superkingdom
2157	Archaea	superkingdom
"""
cumulative_groups = ["2864", "33630", "2830", "33634", "2836", "39119", "33154", "35675", "3041", "33208", "543769", "431838", "91989", "3027", "2", "5878", "10239", "2157"]


# and we'll also need to distinguish a few nested groups.
"""
Alveolata - Dinophyceae, Ciliophora
Stramenopiles - Bacillariophyta, Dictyochophyceae, Pelagophyceae, Bolidophyceae
Opisthokonta - Metazoa
"""
# we'll subtract the counts in the values from the counts in the key (correspond to groups above)
adjusted_groups = { "33630":["2864","5878"], "33634":["2836","39119","35675","91989"], "33154":["33208"] }


# go through and parse the taxa_csv file
TaxaDict = {}
CountsDict = {}

# populate CountsDict with the tax_ids in the cumulative_groups and discrete_groups lists:
for tax_id in discrete_groups:
	CountsDict[tax_id] = {"total_counts":0}
for tax_id in cumulative_groups:
	CountsDict[tax_id] = {"total_counts":0}
# and we need a bin to catch 'Other'
CountsDict["0"] = {"total_counts":0, "rank":"NA", "tax_name":"Other"}

taxa_csv = csv.DictReader(open((args.taxa_csv),'r'))
# each tax_id gets its own dictionary from the row
for row in taxa_csv:
	tax_id = row["tax_id"]
	TaxaDict[tax_id] = row
	# also, add rank and name to the CountsDict:
	if tax_id in CountsDict:
		CountsDict[tax_id]["rank"] = row["rank"]
		CountsDict[tax_id]["tax_name"] = row["tax_name"]

def count_discrete(tax_id, count):

	CountsDict[tax_id]["total_counts"] += count

def count_cumulative(tax_id, count):

	ranks = ["root","below_root","superkingdom","below_superkingdom","below_below_superkingdom","below_below_below_superkingdom","below_below_below_below_superkingdom","kingdom","below_kingdom","below_below_kingdom","below_below_below_kingdom","below_below_below_below_kingdom","below_below_below_below_below_kingdom","subkingdom","phylum","below_phylum","below_below_phylum","below_below_below_phylum","below_below_below_below_phylum","below_below_below_below_below_phylum","below_below_below_below_below_below_phylum","below_below_below_below_below_below_below_phylum","below_below_below_below_below_below_below_below_phylum","below_below_below_below_below_below_below_below_below_phylum","below_below_below_below_below_below_below_below_below_below_phylum","subphylum","below_subphylum","below_below_subphylum","below_below_below_subphylum","below_below_below_below_subphylum","below_below_below_below_below_subphylum","below_below_below_below_below_below_subphylum","below_below_below_below_below_below_below_subphylum","below_below_below_below_below_below_below_below_subphylum","superclass","class","below_class","below_below_class","below_below_below_class","subclass","below_subclass","infraclass","below_infraclass","below_below_infraclass","below_below_below_infraclass","below_below_below_below_infraclass","superorder","below_superorder","order","below_order","below_below_order","suborder","superfamily","family","below_family","below_below_family","subfamily","tribe","genus","below_genus","species_group","species","below_species","below_below_species","subspecies","below_subspecies","varietas","below_varietas","forma"]

	# and we'll add that count to each of the parent ranks (and the identity rank)
	if tax_id in TaxaDict:
		for rank in ranks:
			rank_id = TaxaDict[tax_id][rank]
			# here is where we increment the count:
			if len(rank_id) > 0 and rank_id in cumulative_groups:
				CountsDict[rank_id]["total_counts"] += count


# then, load in the counts file and iterate through the taxids
counts_data = csv.DictReader(open((args.csv),'r'))
# now we'll go through each tax_id counts row:
ranks = ["root","below_root","superkingdom","below_superkingdom","below_below_superkingdom","below_below_below_superkingdom","below_below_below_below_superkingdom","kingdom","below_kingdom","below_below_kingdom","below_below_below_kingdom","below_below_below_below_kingdom","below_below_below_below_below_kingdom","subkingdom","phylum","below_phylum","below_below_phylum","below_below_below_phylum","below_below_below_below_phylum","below_below_below_below_below_phylum","below_below_below_below_below_below_phylum","below_below_below_below_below_below_below_phylum","below_below_below_below_below_below_below_below_phylum","below_below_below_below_below_below_below_below_below_phylum","below_below_below_below_below_below_below_below_below_below_phylum","subphylum","below_subphylum","below_below_subphylum","below_below_below_subphylum","below_below_below_below_subphylum","below_below_below_below_below_subphylum","below_below_below_below_below_below_subphylum","below_below_below_below_below_below_below_subphylum","below_below_below_below_below_below_below_below_subphylum","superclass","class","below_class","below_below_class","below_below_below_class","subclass","below_subclass","infraclass","below_infraclass","below_below_infraclass","below_below_below_infraclass","below_below_below_below_infraclass","superorder","below_superorder","order","below_order","below_below_order","suborder","superfamily","family","below_family","below_below_family","subfamily","tribe","genus","below_genus","species_group","species","below_species","below_below_species","subspecies","below_subspecies","varietas","below_varietas","forma"]
for row in counts_data:
	tax_id = row["tax_id"]
	# for each 'sample_id', we will add the count to the tax_id in the field for any of the ranks:
	count = int(row["count"].strip())
	COUNTS_KEPT = False
	if tax_id in discrete_groups:
		count_discrete(tax_id, count)
	elif tax_id in TaxaDict:
		for rank in ranks:
			rank_id = TaxaDict[tax_id][rank]
			# here is where we increment the count:
			if len(rank_id) > 0 and rank_id in cumulative_groups:
				CountsDict[rank_id]["total_counts"] += count
				COUNTS_KEPT = True
			elif len(rank_id) > 0 and rank_id not in cumulative_groups:
				CountsDict["0"]["total_counts"] += count
	else:
		CountsDict["0"]["total_counts"] += count

# finally, we'll adjust the counts for some of the nested groups we want to keep separate:
for tax_id in adjusted_groups.keys():
	for nested_tax_id in adjusted_groups[tax_id]:
		CountsDict[tax_id]["total_counts"] -= CountsDict[nested_tax_id]["total_counts"]


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
