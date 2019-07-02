#!/usr/bin/env python



# d1.sum_by_group_and_pfam_FPKMx.py
"""
# NEW FOR d1.sum_by_group_and_pfam2.py:
# compute FPKMx-normalized values as well for each time point
# (previous version only does summed standard-normalized counts)

# counts summed up based on pfam_id and phylo_group (e.g., all counts for haptophyte ferritin)
# we need to sum up THREE NUMBERS for EACH of 24 time points:
# max, min, and mean (for our n=2 samples)

# and we can compute more advanced statistics (stdev? sde?) for the 8 points within each time

# input csv structure:
tax_id,pfam_id,phylo_group,num_contigs,S11C1_A_1800,S11C1_C_1800,S14C1_B_2200,S14C1_C_2200,S15C1_B_200,S15C1_C_200,
S16C1_A_600,S16C1_B_600,S17C1_A_1000,S17C1_B_1000,S18C1_A_1400,S18C1_C_1400,S19C1_A_1800,S19C1_C_1800,S20C1_B_2200,
S20C1_C_2200,S21C1_B_200,S21C1_C_200,S22C1_A_600,S22C1_C_600,S23C1_B_1000,S23C1_C_1000,S24C1_A_1400,S24C1_B_1400,
S26C1_A_1800,S26C1_C_1800,S28C1_B_2200,S28C1_C_2200,S29C1_A_200,S29C1_C_200,S30C1_A_600,S30C1_C_600,S31C1_A_1000,
S31C1_C_1000,S32C1_B_1400,S32C1_C_1400,S33C1_A_1800,S33C1_C_1800,S34C1_B_2200,S34C1_C_2200,S35C1_A_200,S35C1_C_200,
S06C1_A_600,S06C1_C_600,S07C1_A_1000,S07C1_B_1000,S08C1_B_1400,S08C1_C_1400

# output csv structure (sum of FPKMx-normalized for all sample ids)
phylo_group,pfam_id,S11C1_A_1800,S11C1_C_1800,S14C1_B_2200,S14C1_C_2200,S15C1_B_200,S15C1_C_200,
S16C1_A_600,S16C1_B_600,S17C1_A_1000,S17C1_B_1000,S18C1_A_1400,S18C1_C_1400,S19C1_A_1800,S19C1_C_1800,S20C1_B_2200,
S20C1_C_2200,S21C1_B_200,S21C1_C_200,S22C1_A_600,S22C1_C_600,S23C1_B_1000,S23C1_C_1000,S24C1_A_1400,S24C1_B_1400,
S26C1_A_1800,S26C1_C_1800,S28C1_B_2200,S28C1_C_2200,S29C1_A_200,S29C1_C_200,S30C1_A_600,S30C1_C_600,S31C1_A_1000,
S31C1_C_1000,S32C1_B_1400,S32C1_C_1400,S33C1_A_1800,S33C1_C_1800,S34C1_B_2200,S34C1_C_2200,S35C1_A_200,S35C1_C_200,
S06C1_A_600,S06C1_C_600,S07C1_A_1000,S07C1_B_1000,S08C1_B_1400,S08C1_C_1400

"""
# total group counts, M:
M_path = "/Users/rgroussman/data/SCOPE/diel1/assemblies/kallisto/FPKMx/summed_counts/d1.cumulative_counts_by_tax_id.csv"

sample_ids = ["S11C1_A_1800","S11C1_C_1800","S14C1_B_2200","S14C1_C_2200","S15C1_B_200","S15C1_C_200","S16C1_A_600","S16C1_B_600","S17C1_A_1000","S17C1_B_1000","S18C1_A_1400","S18C1_C_1400","S19C1_A_1800","S19C1_C_1800","S20C1_B_2200","S20C1_C_2200","S21C1_B_200","S21C1_C_200","S22C1_A_600","S22C1_C_600","S23C1_B_1000","S23C1_C_1000","S24C1_A_1400","S24C1_B_1400","S26C1_A_1800","S26C1_C_1800","S28C1_B_2200","S28C1_C_2200","S29C1_A_200","S29C1_C_200","S30C1_A_600","S30C1_C_600","S31C1_A_1000","S31C1_C_1000","S32C1_B_1400","S32C1_C_1400","S33C1_A_1800","S33C1_C_1800","S34C1_B_2200","S34C1_C_2200","S35C1_A_200","S35C1_C_200","S06C1_A_600","S06C1_C_600","S07C1_A_1000","S07C1_B_1000","S08C1_B_1400","S08C1_C_1400"]

import argparse
import csv
#import numpy as np

# accept arguments:
parser = argparse.ArgumentParser()
parser.add_argument("-c", "--csv", help="Input all_data.csv", type=str)
parser.add_argument("-o", "--out_file", help="Specify the outfile ", type=str)
parser.add_argument("-t", "--taxa_csv", help="Input taxa.csv", type=str)
parser.add_argument("-p", "--phylo_groups", help="CSV file with taxonomies to watch out for", type=str)
args = parser.parse_args()

def init_taxa_dict(tax_data):

	TaxDict = {}
	for row in tax_data:
		tax_id = row['tax_id']
		TaxDict[tax_id] = row
	return TaxDict

def init_phylo_groups(phylo_groups):

	# load the file and populate a dictionary
	# with data for the specific tax groups to look for
	# header: name,tax_id,rank,hi_group
	PhyloDict = {}
	for row in phylo_groups:
		tax_id = row['tax_id']
		PhyloDict[tax_id] = {'name':row['name'], 'tax_id':row['tax_id'], 'rank':row['rank'], 'hi_group':row['hi_group']}
	return PhyloDict

def fix_sample_id(sample_name):

	needs_fixing_dict = {'S6C1_A_600':'S06C1_A_600', 'S6C1_C_600':'S06C1_C_600', 'S7C1_A_1000':'S07C1_A_1000', 'S7C1_B_1000':'S07C1_B_1000', 'S8C1_B_1400':'S08C1_B_1400', 'S8C1_C_1400':'S08C1_C_1400'}
	if sample_name in needs_fixing_dict.keys():
		return needs_fixing_dict[sample_name]
	else:
		return sample_name

def check_major_phylo_groups(TaxDict, tax_id):

	# Here, we want to see if this contig has membership in
	# taxonomic groups of particular interest:

	# these are mutually exclusive groups, so we can check each one then:

	if tax_id in TaxDict.keys():
		for phylo_id in PhyloDict.keys():
			if TaxDict[tax_id][PhyloDict[phylo_id]['rank']] == PhyloDict[phylo_id]['tax_id']:
				return PhyloDict[phylo_id]['name'],PhyloDict[phylo_id]['tax_id'],PhyloDict[phylo_id]['hi_group']
	return "NA","NA","NA"

def parse_all_dat_row(CountsDict, row):

	sample_ids = ["S11C1_A_1800","S11C1_C_1800","S14C1_B_2200","S14C1_C_2200","S15C1_B_200","S15C1_C_200","S16C1_A_600","S16C1_B_600","S17C1_A_1000","S17C1_B_1000","S18C1_A_1400","S18C1_C_1400","S19C1_A_1800","S19C1_C_1800","S20C1_B_2200","S20C1_C_2200","S21C1_B_200","S21C1_C_200","S22C1_A_600","S22C1_C_600","S23C1_B_1000","S23C1_C_1000","S24C1_A_1400","S24C1_B_1400","S26C1_A_1800","S26C1_C_1800","S28C1_B_2200","S28C1_C_2200","S29C1_A_200","S29C1_C_200","S30C1_A_600","S30C1_C_600","S31C1_A_1000","S31C1_C_1000","S32C1_B_1400","S32C1_C_1400","S33C1_A_1800","S33C1_C_1800","S34C1_B_2200","S34C1_C_2200","S35C1_A_200","S35C1_C_200","S06C1_A_600","S06C1_C_600","S07C1_A_1000","S07C1_B_1000","S08C1_B_1400","S08C1_C_1400"]

	pfam_id = row['pfam_id']
	tax_id = row['tax_id']
	phylo_group,phylo_id,hi_group = check_major_phylo_groups(TaxDict, tax_id)
	# now, if phylo_group = NA we can go ahead and skip everything b/c we don't care about these nobodies:
	if phylo_group == "NA":
		return
	# add phylo_group as first level in the dict:
	if phylo_group not in CountsDict.keys():
		CountsDict[phylo_group] = {}
		CountsDict[phylo_group]['hi_group'] = hi_group
	if pfam_id not in CountsDict[phylo_group].keys():
		CountsDict[phylo_group][pfam_id] = {"S11C1_A_1800":0, "S11C1_C_1800":0, "S14C1_B_2200":0, "S14C1_C_2200":0, "S15C1_B_200":0, "S15C1_C_200":0, "S16C1_A_600":0, "S16C1_B_600":0, "S17C1_A_1000":0, "S17C1_B_1000":0, "S18C1_A_1400":0, "S18C1_C_1400":0, "S19C1_A_1800":0, "S19C1_C_1800":0, "S20C1_B_2200":0, "S20C1_C_2200":0, "S21C1_B_200":0, "S21C1_C_200":0, "S22C1_A_600":0, "S22C1_C_600":0, "S23C1_B_1000":0, "S23C1_C_1000":0, "S24C1_A_1400":0, "S24C1_B_1400":0, "S26C1_A_1800":0, "S26C1_C_1800":0, "S28C1_B_2200":0, "S28C1_C_2200":0, "S29C1_A_200":0, "S29C1_C_200":0, "S30C1_A_600":0, "S30C1_C_600":0, "S31C1_A_1000":0, "S31C1_C_1000":0, "S32C1_B_1400":0, "S32C1_C_1400":0, "S33C1_A_1800":0, "S33C1_C_1800":0, "S34C1_B_2200":0, "S34C1_C_2200":0, "S35C1_A_200":0, "S35C1_C_200":0, "S06C1_A_600":0, "S06C1_C_600":0, "S07C1_A_1000":0, "S07C1_B_1000":0, "S08C1_B_1400":0, "S08C1_C_1400":0}
	# now we can start adding in times:
	for sample_name in sample_ids:
		# and then we can grab the value of M (millions of bases) for the relevant phylo_group
		Mx = (float(MDict[phylo_id][sample_name]) / 1e6) # and divide by a million to get M
		# now, instead of just adding the counts, we will compute FpKM first (fragments per kilobase per million reads)
		# FpK has already been computed, so we need only divide the counts by Mx
		if Mx == 0:
			FpKMx = 0
		elif Mx > 0:
			FpKMx = float(row[sample_name]) / Mx
			# and then we can add this to the relevant dict slot:
		CountsDict[phylo_group][pfam_id][sample_name] += FpKMx


def parse_M_data(MDict, row):

	tax_id = row['tax_id']
	MDict[tax_id] = row


def output_csv(CountsDict, out_file_name):

	# now we'll open up our output file and build the header:
	output_csv = open(args.out_file, 'w')
	header = ",".join(["phylo_group", "pfam_id", "hi_group", "S11C1_A_1800","S11C1_C_1800","S14C1_B_2200","S14C1_C_2200","S15C1_B_200","S15C1_C_200","S16C1_A_600","S16C1_B_600","S17C1_A_1000","S17C1_B_1000","S18C1_A_1400","S18C1_C_1400","S19C1_A_1800","S19C1_C_1800","S20C1_B_2200","S20C1_C_2200","S21C1_B_200","S21C1_C_200","S22C1_A_600","S22C1_C_600","S23C1_B_1000","S23C1_C_1000","S24C1_A_1400","S24C1_B_1400","S26C1_A_1800","S26C1_C_1800","S28C1_B_2200","S28C1_C_2200","S29C1_A_200","S29C1_C_200","S30C1_A_600","S30C1_C_600","S31C1_A_1000","S31C1_C_1000","S32C1_B_1400","S32C1_C_1400","S33C1_A_1800","S33C1_C_1800","S34C1_B_2200","S34C1_C_2200","S35C1_A_200","S35C1_C_200","S06C1_A_600","S06C1_C_600","S07C1_A_1000","S07C1_B_1000","S08C1_B_1400","S08C1_C_1400"])
	output_csv.write(header + "\n")

	# and go through the big dictionary and write out the counts:
	for phylo_group in CountsDict.keys():
		hi_group = CountsDict[phylo_group]['hi_group']
		for pfam_id in CountsDict[phylo_group].keys():
			if pfam_id == 'hi_group':
				continue
			outline_elts = [phylo_group, pfam_id, hi_group]
			for sample_name in sample_ids:
				outline_elts.append(CountsDict[phylo_group][pfam_id][sample_name])
			outline = ",".join(str(elt) for elt in outline_elts) + "\n"
			output_csv.write(outline)


# step 0b: Load in 'M' counts file:
print "Loading in Mx counts..."
MDict = {}
M_data = csv.DictReader(open(M_path))
for row in M_data:
	parse_M_data(MDict, row)


# step 0c: load in taxa.csv file and parse
print "Loading and parsing NCBI taxonomy..."
tax_data = csv.DictReader(open(args.taxa_csv))
TaxDict = init_taxa_dict(tax_data)

# step 0d: load in phylo_groups csv
print "Loading and parsing phylo groups..."
phylo_groups = csv.DictReader(open(args.phylo_groups))
PhyloDict = init_phylo_groups(phylo_groups)


# step 1: load in an 'd1.pfam_counts_by_phylo_group.csv' file:
print "Loading all counts..."
CountsDict = {}
all_data = csv.DictReader(open(args.csv))
for row in all_data:
	parse_all_dat_row(CountsDict, row)

# step 3: output results to a new csv file:
print "Writing out results..."
output_csv(CountsDict, args.out_file)
print "Finished!"
