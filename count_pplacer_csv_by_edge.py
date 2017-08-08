#!/usr/bin/env python

# count_pplacer_csv_by_edge.py

"""
Ryan Groussman
Armbrust Lab, University of Washington 2017

This script will take as input the 'guppy to_csv' CSV output from pplacer jplace
files. Going through the list of pqueries (at the moment designed for the 'keep
at most 1' option) and will count the incidence of each EDGE
"""

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-l", "--list_csv", help="Specify path to a list of guppy csv files.", type=str)
parser.add_argument("-n", "--normfactors_csv", help="Specify path to an csv file with norm factors for your samples.", type=str)
parser.add_argument("-o", "--out_file", help="Name of output csv", type=str)
args = parser.parse_args()


def initialize_norm_counts_dict():
	"""from
	default is diel1_standard_counts_SUMS.txt (via Bryn Durham, 13Jan17)
	Norm factor is the sum of BPDSmall1, BPDSmall2, BPDMedium1, BPDMedium2,
	BPDLarge1 and BPDLarge2 from bowtie2 counts, then averaged, and used to
	divide total 'added' counts along with volume:
	NORM FACTOR = ADDED / AVG / VOLUME

	If called with -n or --normfactors_csv, will load a CSV with norm factors
	in this format specified to your samples, with sample names and machine
	reps in the first field separated by a period:

		S14C1_0_2umC.NS7,4779.004092
		S14C1_0_2umC.NS8,4891.099897
		S14C1_3umA.NS1,51153.58824
		S14C1_3umA.NS2,17896.11284


	If called with --machine_runs, machine run data is kept separate.
	Otherwise the machine run data is averaged together.
	"""


	if args.normfactors_csv != None:
		NormFactorsDict = {}


		if args.machine_runs == True:
			# keep machine runs separate
			normfactors_csv = open(args.normfactors_csv, 'r')
			for line in normfactors_csv:
				line_elts = line.split(",")
				run = line_elts[0]
				norm_factor = float(line_elts[1].strip())
				NormFactorsDict[run] = norm_factor
			return NormFactorsDict

		elif args.machine_runs == False:
			# collect them all
			FullNormFactorsDict = {}
			normfactors_csv = open(args.normfactors_csv, 'r')
			for line in normfactors_csv:
				line_elts = line.split(",")
				sample = line_elts[0].split(".")[0]
				machine_run = line_elts[0].split(".")[1]
				norm_factor = float(line_elts[1].strip())
				if sample not in FullNormFactorsDict:
					FullNormFactorsDict[sample] = {}
				FullNormFactorsDict[sample][machine_run] = norm_factor

			# now average the normfactors for each machine_run in a sample:
			for sample in FullNormFactorsDict:
				norm_factor_sum = 0
				mr_count = 0
				for machine_run in FullNormFactorsDict[sample]:
					norm_factor_sum += FullNormFactorsDict[sample][machine_run]
					mr_count += 1
				avg_norm_factor = norm_factor_sum / mr_count
				NormFactorsDict[sample] = avg_norm_factor

			return NormFactorsDict

		return NormFactorsDict
	elif args.normfactors_csv == None:
		# ( to maintain back compatibility - these are the 'default' norm factors if not called with --normfactors_csv )
		# Values taken from diel1_standard_counts_SUMS.txt (via Bryn Durham, 13Jan17)
		NormFactorsDict = {
		'S11C1_A_1800':2393,
		'S11C1_C_1800':2254,
		'S14C1_B_2200':2978,
		'S14C1_C_2200':2365,
		'S15C1_B_200':2674,
		'S15C1_C_200':1970,
		'S16C1_A_600':2597,
		'S16C1_B_600':2520,
		'S17C1_A_1000':2270,
		'S17C1_B_1000':1904,
		'S18C1_A_1400':1339,
		'S18C1_C_1400':2287,
		'S19C1_A_1800':2485,
		'S19C1_C_1800':2097,
		'S20C1_B_2200':2506,
		'S20C1_C_2200':2666,
		'S21C1_B_200':2107,
		'S21C1_C_200':2172,
		'S22C1_A_600':2059,
		'S22C1_C_600':1654,
		'S23C1_B_1000':3005,
		'S23C1_C_1000':1577,
		'S24C1_A_1400':1917,
		'S24C1_B_1400':1466,
		'S26C1_A_1800':2833,
		'S26C1_C_1800':2327,
		'S28C1_B_2200':1247,
		'S28C1_C_2200':2669,
		'S29C1_A_200':2240,
		'S29C1_C_200':1348,
		'S30C1_A_600':1831,
		'S30C1_C_600':2598,
		'S31C1_A_1000':2192,
		'S31C1_C_1000':1834,
		'S32C1_B_1400':1791,
		'S32C1_C_1400':2055,
		'S33C1_A_1800':1949,
		'S33C1_C_1800':2500,
		'S34C1_B_2200':1858,
		'S34C1_C_2200':2641,
		'S35C1_A_200':2625,
		'S35C1_C_200':2630,
		'S6C1_A_600':1396,
		'S6C1_C_600':3435,
		'S7C1_A_1000':1777,
		'S7C1_B_1000':1560,
		'S8C1_B_1400':1887,
		'S8C1_C_1400':4185
		}

	return NormFactorsDict

def process_guppy_csv(csv_file):
	"""Process a guppy_csv"""

	total_pqueries = 0
	sample_name = csv_file.split(".")[1]

	current_csv = open(csv_file, 'r')
	current_csv.readline()
	for line in current_csv:
		total_pqueries += 1
		line_elts = line.split(",")
		edge_num = line_elts[3]
		if edge_num not in CountsDict:
			CountsDict[edge_num] = {}
		if sample_name not in CountsDict[edge_num]:
			CountsDict[edge_num][sample_name] = 1
		elif sample_name in CountsDict[edge_num]:
			CountsDict[edge_num][sample_name] += 1
	CountsDict["Total"][sample_name] = total_pqueries
	current_csv.close()


csv_list = open(args.list_csv, 'r')
out_file = open(args.out_file, 'w')

NormFactorsDict = initialize_norm_counts_dict()
# print NormFactorsDict

CountsDict = {}
CountsDict["Total"] = {}

sample_list = []
for csv_file in csv_list:
	csv_file = csv_file.strip()
	process_guppy_csv(csv_file)
	sample_name = csv_file.split(".")[1]
	sample_list.append(sample_name)


# after counting, multiply by norm factor and write out results:
header = "edge_num," + ",".join(sample_list)  + "\n"
out_file.write(header)
for edge in CountsDict:
	edge_norm_counts_list = []
	for sample_name in sample_list:
		if sample_name in CountsDict[edge]:
			norm_factor = NormFactorsDict[sample_name]
			norm_count = CountsDict[edge][sample_name] * norm_factor
		else:
			norm_count = 0
		edge_norm_counts_list.append(str(norm_count))
	out_string = str(edge) + "," + ",".join(edge_norm_counts_list) + "\n"
	out_file.write(out_string)

# write out totals:
# totals_list = ["Total"]
# for sample_name in sample_list:
# 	norm_factor = NormFactorsDict[sample_name]
# 	totals_list.append(str(CountsDict["Total"][sample_name] * norm_factor))
# out_string = ",".join(totals_list) + "\n"
# out_file.write(out_string)
