#!/usr/bin/env python

# d1.cumulative_count_sums.py

"""
This script will calculate the cumulative, nested sums for the
combined counts file output by d1.combine_count_sums.py

We will leverage nested taxonomic information to determine the cumulative sums
for each tax_id inclusive of its children tax_ids

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
	CountsDict[tax_id] = {"total_counts":0,"S11C1_A_1800":0, "S11C1_C_1800":0, "S14C1_B_2200":0, "S14C1_C_2200":0, "S15C1_B_200":0, "S15C1_C_200":0, "S16C1_A_600":0, "S16C1_B_600":0, "S17C1_A_1000":0, "S17C1_B_1000":0, "S18C1_A_1400":0, "S18C1_C_1400":0, "S19C1_A_1800":0, "S19C1_C_1800":0, "S20C1_B_2200":0, "S20C1_C_2200":0, "S21C1_B_200":0, "S21C1_C_200":0, "S22C1_A_600":0, "S22C1_C_600":0, "S23C1_B_1000":0, "S23C1_C_1000":0, "S24C1_A_1400":0, "S24C1_B_1400":0, "S26C1_A_1800":0, "S26C1_C_1800":0, "S28C1_B_2200":0, "S28C1_C_2200":0, "S29C1_A_200":0, "S29C1_C_200":0, "S30C1_A_600":0, "S30C1_C_600":0, "S31C1_A_1000":0, "S31C1_C_1000":0, "S32C1_B_1400":0, "S32C1_C_1400":0, "S33C1_A_1800":0, "S33C1_C_1800":0, "S34C1_B_2200":0, "S34C1_C_2200":0, "S35C1_A_200":0, "S35C1_C_200":0, "S06C1_A_600":0, "S06C1_C_600":0, "S07C1_A_1000":0, "S07C1_B_1000":0, "S08C1_B_1400":0, "S08C1_C_1400":0}
	CountsDict[tax_id]["rank"] = row["rank"]
	CountsDict[tax_id]["tax_name"] = row["tax_name"]


# create the CountsDict

# then, load in the counts file and iterate through the taxids
counts_data = csv.DictReader(open((args.csv),'r'))
sample_ids = ["total_counts","S11C1_A_1800","S11C1_C_1800","S14C1_B_2200","S14C1_C_2200","S15C1_B_200","S15C1_C_200","S16C1_A_600","S16C1_B_600","S17C1_A_1000","S17C1_B_1000","S18C1_A_1400","S18C1_C_1400","S19C1_A_1800","S19C1_C_1800","S20C1_B_2200","S20C1_C_2200","S21C1_B_200","S21C1_C_200","S22C1_A_600","S22C1_C_600","S23C1_B_1000","S23C1_C_1000","S24C1_A_1400","S24C1_B_1400","S26C1_A_1800","S26C1_C_1800","S28C1_B_2200","S28C1_C_2200","S29C1_A_200","S29C1_C_200","S30C1_A_600","S30C1_C_600","S31C1_A_1000","S31C1_C_1000","S32C1_B_1400","S32C1_C_1400","S33C1_A_1800","S33C1_C_1800","S34C1_B_2200","S34C1_C_2200","S35C1_A_200","S35C1_C_200","S06C1_A_600","S06C1_C_600","S07C1_A_1000","S07C1_B_1000","S08C1_B_1400","S08C1_C_1400"]
ranks = ["root","below_root","superkingdom","below_superkingdom","below_below_superkingdom","below_below_below_superkingdom","below_below_below_below_superkingdom","kingdom","below_kingdom","below_below_kingdom","below_below_below_kingdom","below_below_below_below_kingdom","below_below_below_below_below_kingdom","subkingdom","phylum","below_phylum","below_below_phylum","below_below_below_phylum","below_below_below_below_phylum","below_below_below_below_below_phylum","below_below_below_below_below_below_phylum","below_below_below_below_below_below_below_phylum","below_below_below_below_below_below_below_below_phylum","below_below_below_below_below_below_below_below_below_phylum","below_below_below_below_below_below_below_below_below_below_phylum","subphylum","below_subphylum","below_below_subphylum","below_below_below_subphylum","below_below_below_below_subphylum","below_below_below_below_below_subphylum","below_below_below_below_below_below_subphylum","below_below_below_below_below_below_below_subphylum","below_below_below_below_below_below_below_below_subphylum","superclass","class","below_class","below_below_class","below_below_below_class","subclass","below_subclass","infraclass","below_infraclass","below_below_infraclass","below_below_below_infraclass","below_below_below_below_infraclass","superorder","below_superorder","order","below_order","below_below_order","suborder","superfamily","family","below_family","below_below_family","subfamily","tribe","genus","below_genus","species_group","species","below_species","below_below_species","subspecies","below_subspecies","varietas","below_varietas","forma"]
# now we'll go through each tax_id counts row:
for row in counts_data:
	tax_id = row["tax_id"]
	# for each 'sample_id', we will add the count to the tax_id in the field for any of the ranks:
	for sample in sample_ids:
		count = float(row[sample].strip())
		# and we'll add that count to each of the parent ranks (and the identity rank)
		for rank in ranks:
			if tax_id in TaxaDict:
				rank_id = TaxaDict[tax_id][rank]
			else:
				rank_id = ""
			# here is where we increment the count:
			if len(rank_id) > 0:
				CountsDict[rank_id][sample] += count

# output:
output_csv = open(args.output, 'w')
# this is the same ordering as sample_ids above, but with tax_id, and S6C1 -> S06C1 fixed, etc
header = ",".join(["tax_id","rank","tax_name","total_counts","S11C1_A_1800","S11C1_C_1800","S14C1_B_2200","S14C1_C_2200","S15C1_B_200","S15C1_C_200","S16C1_A_600","S16C1_B_600","S17C1_A_1000","S17C1_B_1000","S18C1_A_1400","S18C1_C_1400","S19C1_A_1800","S19C1_C_1800","S20C1_B_2200","S20C1_C_2200","S21C1_B_200","S21C1_C_200","S22C1_A_600","S22C1_C_600","S23C1_B_1000","S23C1_C_1000","S24C1_A_1400","S24C1_B_1400","S26C1_A_1800","S26C1_C_1800","S28C1_B_2200","S28C1_C_2200","S29C1_A_200","S29C1_C_200","S30C1_A_600","S30C1_C_600","S31C1_A_1000","S31C1_C_1000","S32C1_B_1400","S32C1_C_1400","S33C1_A_1800","S33C1_C_1800","S34C1_B_2200","S34C1_C_2200","S35C1_A_200","S35C1_C_200","S06C1_A_600","S06C1_C_600","S07C1_A_1000","S07C1_B_1000","S08C1_B_1400","S08C1_C_1400"])
output_csv.write(header + "\n")

for tax_id in CountsDict.keys():
	outline_elts = [str(tax_id)]
	outline_elts.append(str(CountsDict[tax_id]["rank"]))
	outline_elts.append(str(CountsDict[tax_id]["tax_name"]))
	for sample in sample_ids:
		outline_elts.append(str(CountsDict[tax_id][sample]))
	output_csv.write(",".join(outline_elts) + "\n")



# for each sample (and total counts) for each
# add to the cumulative counts for each tax_id found in the ranks to which it belongs.
# e.g.: for tax_id 33630 'Alveolata', add the counts to below_root 131567	superkingdom 2759	below_superkingdom 33630
# so counts are retained at its own discrete level
