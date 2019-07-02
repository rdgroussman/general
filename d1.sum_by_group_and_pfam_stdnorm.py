#!/usr/bin/env python



# d1.sum_by_group_and_pfam_FPKMx.py
"""

# input:

# output: summed standard-normalized counts for each time point

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
"""

# hard-coded paths:
time_dict_path = "/Users/rgroussman/data/SCOPE/diel1/diel1_timepoints_wreps.csv"

# total group counts, M (for FPKMx, not necessary here)
#M_path = "/Users/rgroussman/data/SCOPE/diel1/assemblies/kallisto/FPKMx/summed_counts/d1.cumulative_counts_by_tax_id.csv"


import argparse
import csv
#import numpy as np

# accept arguments:
parser = argparse.ArgumentParser()
parser.add_argument("-c", "--csv", help="Input all_data.csv", type=str)
parser.add_argument("-o", "--out_file_prefix", help="Specify the outfile prefix", type=str)
parser.add_argument("-p", "--phylo_groups", help="CSV file with taxonomies to watch out for", type=str)
parser.add_argument("-t", "--taxa_csv", help="Input taxa.csv", type=str)
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

def get_times(sample_name):

	return TimeDict[sample_name]["clock_time"],TimeDict[sample_name]["abs_time"],TimeDict[sample_name]["rep"]

def parse_all_dat_row(CountsDict, row):

	samples = ["S11C1_A_1800","S11C1_C_1800","S14C1_B_2200","S14C1_C_2200","S15C1_B_200","S15C1_C_200","S16C1_A_600","S16C1_B_600","S17C1_A_1000","S17C1_B_1000","S18C1_A_1400","S18C1_C_1400","S19C1_A_1800","S19C1_C_1800","S20C1_B_2200","S20C1_C_2200","S21C1_B_200","S21C1_C_200","S22C1_A_600","S22C1_C_600","S23C1_B_1000","S23C1_C_1000","S24C1_A_1400","S24C1_B_1400","S26C1_A_1800","S26C1_C_1800","S28C1_B_2200","S28C1_C_2200","S29C1_A_200","S29C1_C_200","S30C1_A_600","S30C1_C_600","S31C1_A_1000","S31C1_C_1000","S32C1_B_1400","S32C1_C_1400","S33C1_A_1800","S33C1_C_1800","S34C1_B_2200","S34C1_C_2200","S35C1_A_200","S35C1_C_200","S06C1_A_600","S06C1_C_600","S07C1_A_1000","S07C1_B_1000","S08C1_B_1400","S08C1_C_1400"]

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
		#CountsDict[phylo_group]['phylo_id'] = phylo_id
	if pfam_id not in CountsDict[phylo_group].keys():
		# assign to rep within a time within a time category (clock or abs)
		CountsDict[phylo_group][pfam_id] = {'clock_time':{2:{1:[],2:[]}, 6:{1:[],2:[]}, 10:{1:[],2:[]}, 14:{1:[],2:[]}, 18:{1:[],2:[]}, 22:{1:[],2:[]}},'abs_time':{0:{1:[],2:[]}, 4:{1:[],2:[]}, 8:{1:[],2:[]}, 12:{1:[],2:[]}, 16:{1:[],2:[]}, 20:{1:[],2:[]}, 24:{1:[],2:[]}, 28:{1:[],2:[]}, 32:{1:[],2:[]}, 36:{1:[],2:[]}, 40:{1:[],2:[]}, 44:{1:[],2:[]}, 48:{1:[],2:[]}, 52:{1:[],2:[]}, 56:{1:[],2:[]}, 60:{1:[],2:[]}, 64:{1:[],2:[]}, 68:{1:[],2:[]}, 72:{1:[],2:[]}, 76:{1:[],2:[]}, 80:{1:[],2:[]}, 84:{1:[],2:[]}, 88:{1:[],2:[]}, 92:{1:[],2:[]}, 92:{1:[],2:[]}, }}
	# now we can start adding in times:
	for sample_name in samples:
		fixed_name = fix_sample_id(sample_name)
		clock_time,abs_time,rep = get_times(fixed_name)
		# now just add the counts to the relevant dict slot:
		CountsDict[phylo_group][pfam_id]['clock_time'][clock_time]
		CountsDict[phylo_group][pfam_id]['clock_time'][clock_time][rep].append(int(row[sample_name]))
		CountsDict[phylo_group][pfam_id]['abs_time'][abs_time][rep].append(int(row[sample_name]))


def compute_stats(CountsDict):

	# now we'll run through the CountsDict and compute stats for each of the group/pfam collections:
	for phylo_group in CountsDict.keys():
		StatsDict[phylo_group] = {}
		for pfam_id in CountsDict[phylo_group].keys():
			if pfam_id == 'hi_group':
				continue
			StatsDict[phylo_group][pfam_id] = {'clock_time':{2:{}, 6:{}, 10:{}, 14:{}, 18:{}, 22:{}},'abs_time':{0:{}, 4:{}, 8:{}, 12:{}, 16:{}, 20:{}, 24:{}, 28:{}, 32:{}, 36:{}, 40:{}, 44:{}, 48:{}, 52:{}, 56:{}, 60:{}, 64:{}, 68:{}, 72:{}, 76:{}, 80:{}, 84:{}, 88:{}, 92:{}, 92:{}, }}
			# clock time:
			for clock_time in CountsDict[phylo_group][pfam_id]['clock_time']:
				#print clock_time
				StatsDict[phylo_group][pfam_id]['clock_time'][clock_time]['sum'] = sum(CountsDict[phylo_group][pfam_id]['clock_time'][clock_time][1]) + sum(CountsDict[phylo_group][pfam_id]['clock_time'][clock_time][2])
				StatsDict[phylo_group][pfam_id]['clock_time'][clock_time]['max'] = max(sum(CountsDict[phylo_group][pfam_id]['clock_time'][clock_time][1]),sum(CountsDict[phylo_group][pfam_id]['clock_time'][clock_time][2]))
				StatsDict[phylo_group][pfam_id]['clock_time'][clock_time]['min'] = min(sum(CountsDict[phylo_group][pfam_id]['clock_time'][clock_time][1]),sum(CountsDict[phylo_group][pfam_id]['clock_time'][clock_time][2]))
				StatsDict[phylo_group][pfam_id]['clock_time'][clock_time]['mean'] = (StatsDict[phylo_group][pfam_id]['clock_time'][clock_time]['sum'] / 2)

			for abs_time in CountsDict[phylo_group][pfam_id]['abs_time']:
				StatsDict[phylo_group][pfam_id]['abs_time'][abs_time]['sum'] = sum(CountsDict[phylo_group][pfam_id]['abs_time'][abs_time][1] + CountsDict[phylo_group][pfam_id]['abs_time'][abs_time][2])
				StatsDict[phylo_group][pfam_id]['abs_time'][abs_time]['max'] = max(sum(CountsDict[phylo_group][pfam_id]['abs_time'][abs_time][1]),sum(CountsDict[phylo_group][pfam_id]['abs_time'][abs_time][2]))
				StatsDict[phylo_group][pfam_id]['abs_time'][abs_time]['min'] = min(sum(CountsDict[phylo_group][pfam_id]['abs_time'][abs_time][1]),sum(CountsDict[phylo_group][pfam_id]['abs_time'][abs_time][2]))
				StatsDict[phylo_group][pfam_id]['abs_time'][abs_time]['mean'] = (StatsDict[phylo_group][pfam_id]['abs_time'][abs_time]['sum'] / 2)


def output_csv(ContigData, out_file_name):

	clock_time_csv = open((out_file_name + ".clock_time.std_norm_data.csv"), 'w')
	header = ['phylo_group', 'pfam_id', 'hi_group', 'time', 'mean', 'min','max','sum']
	clock_time_csv.write(",".join(header) + "\n")
	# first, print out clock_time counts to its own file:
	for phylo_group in CountsDict.keys():
		#print phylo_group
		#print CountsDict[phylo_group]
		hi_group = CountsDict[phylo_group]['hi_group']
		for pfam_id in CountsDict[phylo_group]:
			if pfam_id == 'hi_group':
				continue
			#print pfam_id
			for clock_time in CountsDict[phylo_group][pfam_id]['clock_time']:
				out_fields = [phylo_group, pfam_id, hi_group, str(clock_time)]
				for stat in ['mean', 'min','max','sum']:
					out_fields.append(str(StatsDict[phylo_group][pfam_id]['clock_time'][clock_time][stat]))
				#print out_fields
				#print ",".join(out_fields)
				clock_time_csv.write(",".join(out_fields)+ "\n")

	# and for absolute time:
	abs_time_csv = open((out_file_name + ".abs_time.std_norm_data.csv"), 'w')
	header = ['phylo_group', 'pfam_id', 'hi_group', 'time', 'mean', 'min', 'max', 'sum']
	abs_time_csv.write(",".join(header) + "\n")
	# then, print out abs_time counts:
	for phylo_group in CountsDict.keys():
		hi_group = CountsDict[phylo_group]['hi_group']
		for pfam_id in CountsDict[phylo_group]:
			if pfam_id == 'hi_group':
				continue
			for abs_time in CountsDict[phylo_group][pfam_id]['abs_time']:
				out_fields = [phylo_group, pfam_id, hi_group, str(abs_time)]
				for stat in ['mean', 'min', 'max', 'sum']:
					out_fields.append(str(StatsDict[phylo_group][pfam_id]['abs_time'][abs_time][stat]))
				abs_time_csv.write(",".join(out_fields)+ "\n")

	# and finally, we can also output a FPKMx file:

# step 0a: initialize the TimeDict:
TimeDict = {}
time_csv = csv.DictReader(open(time_dict_path))
for row in time_csv:
	sample_name = row["sample_name"]
	TimeDict[sample_name] = {}
	TimeDict[sample_name]["clock_time"] = int(row["clock_time"])
	TimeDict[sample_name]["abs_time"] = int(row["abs_time"])
	TimeDict[sample_name]["rep"] = int(row["rep"])


# step 0c: load in taxa.csv file and parse
print "Loading and parsing NCBI taxonomy..."
tax_data = csv.DictReader(open(args.taxa_csv))
TaxDict = init_taxa_dict(tax_data)

# step 0d: load in phylo_groups csv
print "Loading and parsing NCBI taxonomy..."
phylo_groups = csv.DictReader(open(args.phylo_groups))
PhyloDict = init_phylo_groups(phylo_groups)

# step 1: load in an 'd1.pfam_counts_by_phylo_group.csv' file:
print "Loading all counts..."
CountsDict = {}
all_data = csv.DictReader(open(args.csv))
for row in all_data:
	parse_all_dat_row(CountsDict, row)

# step 2: compute stats for clock time and abs time
print "Computing stats..."
StatsDict = {}
compute_stats(CountsDict)
#print CountsDict

# step 3: output results to a new csv file:
print "Writing out results..."
output_csv(StatsDict, args.out_file_prefix)
print "Finished!"
