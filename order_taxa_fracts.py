#!/usr/bin/env python

# order_taxa_fracts.py

# This script takes in a taxified taxa_fracts file and will re-organize the output to something like this:
# outputs a CSV file.

"""

name			tax_id	group	S6C1_A_600 S11C1_A_1800 etc...
Prymnesiaceae	234234	hapto	.12			.12			.22
Polarella		234		dino	.01			.01			.33


"""

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="Input taxified_taxa_fracts.csv", type=str)
parser.add_argument("-o", "--out_file", help="Specify the name of the outfile (csv)", type=str)
args = parser.parse_args()

out_csv = open(args.out_file,'w')

def make_time_series_dict():
	"""Links diel1 sample ID to sampling hour"""

	time_series_dict = {
	'S06C1_A_600': "0",
	'S06C1_C_600': "0",
	'S07C1_A_1000': "4",
	'S07C1_B_1000': "4",
	'S08C1_B_1400': "8",
	'S08C1_C_1400': "8",
	'S11C1_A_1800': "12",
	'S11C1_C_1800': "12",
	'S14C1_B_2200': "16",
	'S14C1_C_2200': "16",
	'S15C1_B_200': "20",
	'S15C1_C_200': "20",
	'S16C1_A_600': "24",
	'S16C1_B_600': "24",
	'S17C1_A_1000': "28",
	'S17C1_B_1000': "28",
	'S18C1_A_1400': "32",
	'S18C1_C_1400': "32",
	'S19C1_A_1800': "36",
	'S19C1_C_1800': "36",
	'S20C1_B_2200': "40",
	'S20C1_C_2200': "40",
	'S21C1_B_200': "44",
	'S21C1_C_200': "44",
	'S22C1_A_600': "48",
	'S22C1_C_600': "48",
	'S23C1_B_1000': "52",
	'S23C1_C_1000': "52",
	'S24C1_A_1400': "56",
	'S24C1_B_1400': "56",
	'S26C1_A_1800': "60",
	'S26C1_C_1800': "60",
	'S28C1_B_2200': "64",
	'S28C1_C_2200': "64",
	'S29C1_A_200': "68",
	'S29C1_C_200': "68",
	'S30C1_A_600': "72",
	'S30C1_C_600': "72",
	'S31C1_A_1000': "76",
	'S31C1_C_1000': "76",
	'S32C1_B_1400': "80",
	'S32C1_C_1400': "80",
	'S33C1_A_1800': "84",
	'S33C1_C_1800': "84",
	'S34C1_B_2200': "88",
	'S34C1_C_2200': "88",
	'S35C1_A_200': "92",
	'S35C1_C_200': "92" }

	return time_series_dict

time_series_dict = make_time_series_dict()
TaxaFractsDict = {}

# We'll run through the taxa_fracts.csv and collect the info for each tax_id at each time point.
with open(args.input, 'r') as taxa_csv:
	# next(taxa_csv) # we won't skip the first line (no header)
	for line in taxa_csv:
		line_elts = line.split(",")
		sample = line_elts[0]
		sample_hour = time_series_dict[sample]
		group = line_elts[1]
		tax_id = line_elts[2]
		fraction = line_elts[3]
		name = line_elts[4]
		level = line_elts[5].strip()
		if tax_id not in TaxaFractsDict:
			TaxaFractsDict[tax_id] = {}
			TaxaFractsDict[tax_id]["group"] = group
			TaxaFractsDict[tax_id]["name"] = name
			TaxaFractsDict[tax_id]["level"] = level
		if sample_hour not in TaxaFractsDict[tax_id]:
			TaxaFractsDict[tax_id][sample_hour] = []
		TaxaFractsDict[tax_id][sample_hour].append(fraction)

# Now we'll print out to the out_csv file the info for each tax_id:

ordered_time_points = sorted(set([time_series_dict[key] for key in time_series_dict]), key=int)
header = "tax_id,name,level,group," + ",".join(ordered_time_points) + "\n"

out_csv.write(header)

# print TaxaFractsDict
for tax_id in TaxaFractsDict:
	out_list = [tax_id, TaxaFractsDict[tax_id]["name"], TaxaFractsDict[tax_id]["level"], TaxaFractsDict[tax_id]["group"]]
	for sample_hour in ordered_time_points:
		if sample_hour in TaxaFractsDict[tax_id]:
			avg_fract = sum([float(fract) for fract in TaxaFractsDict[tax_id][sample_hour]]) / float(len(TaxaFractsDict[tax_id][sample_hour]))
		elif sample_hour not in TaxaFractsDict[tax_id]:
			avg_fract = 0
		out_list.append(str(avg_fract))
	out_line = ",".join(out_list) + "\n"
	out_csv.write(out_line)

out_csv.close()
