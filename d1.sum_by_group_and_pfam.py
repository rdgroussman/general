#!/usr/bin/env python



# d1.sum_by_group_and_pfam.py
"""
# what we want:
# counts summed up based on pfam_name and phylo_group (e.g., all counts for haptophyte ferritin)
# we need to sum up THREE NUMBERS for EACH of 24 time points:
# max, min, and mean (for our n=2 samples)

# and we can compute more advanced statistics (stdev? sde?) for the 8 points within each time
"""

# hard-coded paths:
time_dict_path = "/Users/rgroussman/data/SCOPE/diel1/diel1_timepoints.csv"



import argparse
import csv
import numpy as np

# accept arguments:
parser = argparse.ArgumentParser()
parser.add_argument("-c", "--csv", help="Input all_data.csv", type=str)
parser.add_argument("-o", "--out_file_prefix", help="Specify the outfile prefix", type=str)

args = parser.parse_args()


def fix_sample_id(sample_name):

	needs_fixing_dict = {'S6C1_A_600':'S06C1_A_600', 'S6C1_C_600':'S06C1_C_600', 'S7C1_A_1000':'S07C1_A_1000', 'S7C1_B_1000':'S07C1_B_1000', 'S8C1_B_1400':'S08C1_B_1400', 'S8C1_C_1400':'S08C1_C_1400'}
	if sample_name in needs_fixing_dict.keys():
		return needs_fixing_dict[sample_name]
	else:
		return sample_name


def get_times(sample_name):

	return TimeDict[sample_name]["clock_time"],TimeDict[sample_name]["abs_time"]

def parse_all_dat_row(CountsDict, row):

	samples = ["S11C1_A_1800","S11C1_C_1800","S14C1_B_2200","S14C1_C_2200","S15C1_B_200","S15C1_C_200","S16C1_A_600","S16C1_B_600","S17C1_A_1000","S17C1_B_1000","S18C1_A_1400","S18C1_C_1400","S19C1_A_1800","S19C1_C_1800","S20C1_B_2200","S20C1_C_2200","S21C1_B_200","S21C1_C_200","S22C1_A_600","S22C1_C_600","S23C1_B_1000","S23C1_C_1000","S24C1_A_1400","S24C1_B_1400","S26C1_A_1800","S26C1_C_1800","S28C1_B_2200","S28C1_C_2200","S29C1_A_200","S29C1_C_200","S30C1_A_600","S30C1_C_600","S31C1_A_1000","S31C1_C_1000","S32C1_B_1400","S32C1_C_1400","S33C1_A_1800","S33C1_C_1800","S34C1_B_2200","S34C1_C_2200","S35C1_A_200","S35C1_C_200","S06C1_A_600","S06C1_C_600","S07C1_A_1000","S07C1_B_1000","S08C1_B_1400","S08C1_C_1400"]

	nt_id = row['contig_id']
	pfam_name = row['pfam_name']
	phylo_group = row['phylo_group']
	# add phylo_group as first level in the dict:
	if phylo_group not in CountsDict.keys():
		CountsDict[phylo_group] = {}
	if pfam_name not in CountsDict[phylo_group].keys():
		CountsDict[phylo_group][pfam_name] = {'clock_time':{2:[], 6:[], 10:[], 14:[], 18:[], 22:[]},'abs_time':{0:[], 4:[], 8:[], 12:[], 16:[], 20:[], 24:[], 28:[], 32:[], 36:[], 40:[], 44:[], 48:[], 52:[], 56:[], 60:[], 64:[], 68:[], 72:[], 76:[], 80:[], 84:[], 88:[], 92:[], 92:[], }}
	# now we can start adding in times:
	for sample_name in samples:
		fixed_name = fix_sample_id(sample_name)
		clock_time,abs_time = get_times(fixed_name)
		CountsDict[phylo_group][pfam_name]['clock_time'][clock_time].append(int(row[sample_name]))
		CountsDict[phylo_group][pfam_name]['abs_time'][abs_time].append(int(row[sample_name]))

		# print sample_name, clock_time, abs_time
		# CountsDict[phylo_group][pfam_name]['clock']


def compute_stats(CountsDict):

	# now we'll run through the CountsDict and compute stats for each of the group/pfam collections:
	for phylo_group in CountsDict.keys():
		StatsDict[phylo_group] = {}
		for pfam_name in CountsDict[phylo_group]:
			StatsDict[phylo_group][pfam_name] = {'clock_time':{2:{}, 6:{}, 10:{}, 14:{}, 18:{}, 22:{}},'abs_time':{0:{}, 4:{}, 8:{}, 12:{}, 16:{}, 20:{}, 24:{}, 28:{}, 32:{}, 36:{}, 40:{}, 44:{}, 48:{}, 52:{}, 56:{}, 60:{}, 64:{}, 68:{}, 72:{}, 76:{}, 80:{}, 84:{}, 88:{}, 92:{}, 92:{}, }}
			# clock time:
			for clock_time in CountsDict[phylo_group][pfam_name]['clock_time']:
				StatsDict[phylo_group][pfam_name]['clock_time'][clock_time]['sum'] = sum(CountsDict[phylo_group][pfam_name]['clock_time'][clock_time])
				StatsDict[phylo_group][pfam_name]['clock_time'][clock_time]['max'] = max(CountsDict[phylo_group][pfam_name]['clock_time'][clock_time])
				StatsDict[phylo_group][pfam_name]['clock_time'][clock_time]['min'] = min(CountsDict[phylo_group][pfam_name]['clock_time'][clock_time])
				StatsDict[phylo_group][pfam_name]['clock_time'][clock_time]['mean'] = np.mean(CountsDict[phylo_group][pfam_name]['clock_time'][clock_time])
				StatsDict[phylo_group][pfam_name]['clock_time'][clock_time]['stdev'] = np.std(CountsDict[phylo_group][pfam_name]['clock_time'][clock_time])
			for abs_time in CountsDict[phylo_group][pfam_name]['abs_time']:
				StatsDict[phylo_group][pfam_name]['abs_time'][abs_time]['max'] = max(CountsDict[phylo_group][pfam_name]['abs_time'][abs_time])
				StatsDict[phylo_group][pfam_name]['abs_time'][abs_time]['min'] = min(CountsDict[phylo_group][pfam_name]['abs_time'][abs_time])
				StatsDict[phylo_group][pfam_name]['abs_time'][abs_time]['mean'] = np.mean(CountsDict[phylo_group][pfam_name]['abs_time'][abs_time])
				StatsDict[phylo_group][pfam_name]['abs_time'][abs_time]['sum'] = sum(CountsDict[phylo_group][pfam_name]['abs_time'][abs_time])

def output_csv(ContigData, out_file_name):


	clock_time_csv = open((out_file_name + ".clock_time.csv"), 'w')
	header = ['phylo_group', 'pfam_name', 'time', 'mean', 'stdev','sum']
	clock_time_csv.write(",".join(header) + "\n")

	# first, print out clock_time counts to its own file:
	for phylo_group in CountsDict.keys():
		for pfam_name in CountsDict[phylo_group]:
			for clock_time in CountsDict[phylo_group][pfam_name]['clock_time']:
				out_fields = [phylo_group, pfam_name, str(clock_time)]
				for stat in ['mean', 'stdev', 'sum']:
					out_fields.append(str(StatsDict[phylo_group][pfam_name]['clock_time'][clock_time][stat]))
				clock_time_csv.write(",".join(out_fields)+ "\n")

	abs_time_csv = open((out_file_name + ".abs_time.csv"), 'w')
	header = ['phylo_group', 'pfam_name', 'time', 'mean', 'max', 'min', 'sum']
	abs_time_csv.write(",".join(header) + "\n")

	# then, print out abs_time counts:
	for phylo_group in CountsDict.keys():
		for pfam_name in CountsDict[phylo_group]:
			for abs_time in CountsDict[phylo_group][pfam_name]['abs_time']:
				out_fields = [phylo_group, pfam_name, str(abs_time)]
				for stat in ['mean', 'max', 'min', 'sum']:
					out_fields.append(str(StatsDict[phylo_group][pfam_name]['abs_time'][abs_time][stat]))
				abs_time_csv.write(",".join(out_fields)+ "\n")


# step 0: initialize the TimeDict:
TimeDict = {}
with open(time_dict_path, 'r') as time_csv:
	header = next(time_csv) # skip the header:
	for line in time_csv:
		line_elts = line.split(",")
		sample_name = line_elts[0]
		clock_time = int(line_elts[1])
		abs_time = int(line_elts[2].strip())
		TimeDict[sample_name] = {}
		TimeDict[sample_name]["clock_time"] = clock_time
		TimeDict[sample_name]["abs_time"] = abs_time

# step 1: load in an *all_data.csv file:
CountsDict = {}
all_data = csv.DictReader(open(args.csv))
for row in all_data:
	parse_all_dat_row(CountsDict, row)

# step 2: compute stats for clock time and abs time:
StatsDict = {}
compute_stats(CountsDict)
# print CountsDict

# step 3: output results to a new csv file:
output_csv(StatsDict, args.out_file_prefix)
