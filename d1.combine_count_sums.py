#!/usr/bin/env python

# d1.combine_count_sums.py

"""
Refer to notes in calculating_FPKMx.log.sh

This takes as input a list of files, e.g.
	d1.S06C1.est_counts_by_tax_id.csv
	d1.S07C1.est_counts_by_tax_id.csv
	d1.S08C1.est_counts_by_tax_id.csv
	d1.S11C1.est_counts_by_tax_id.csv
	d1.S14C1.est_counts_by_tax_id.csv

These filers were generated by d1.sum_raw_counts.py

These files have the format:
	"tax_id","S11C1_A_1800","S11C1_C_1800","S14C1_B_2200","S14C1_C_2200","S15C1_B_200","S15C1_C_200","S16C1_A_600","S16C1_B_600","S17C1_A_1000","S17C1_B_1000","S18C1_A_1400","S18C1_C_1400","S19C1_A_1800","S19C1_C_1800","S20C1_B_2200","S20C1_C_2200","S21C1_B_200","S21C1_C_200","S22C1_A_600","S22C1_C_600","S23C1_B_1000","S23C1_C_1000","S24C1_A_1400","S24C1_B_1400","S26C1_A_1800","S26C1_C_1800","S28C1_B_2200","S28C1_C_2200","S29C1_A_200","S29C1_C_200","S30C1_A_600","S30C1_C_600","S31C1_A_1000","S31C1_C_1000","S32C1_B_1400","S32C1_C_1400","S33C1_A_1800","S33C1_C_1800","S34C1_B_2200","S34C1_C_2200","S35C1_A_200","S35C1_C_200","S06C1_A_600","S06C1_C_600","S07C1_A_1000","S07C1_B_1000","S08C1_B_1400","S08C1_C_1400"

Output a file with similar format, but with addition of a 'total_counts' field

"""
import argparse
import csv

parser = argparse.ArgumentParser()
parser.add_argument("-l", "--list", help="List of *est_counts_by_tax_id.csv files", type=str)
parser.add_argument("-o", "--output", help="Output CSV path", type=str)
args = parser.parse_args()

def add_file_counts(csv_file):

	counts_data = csv.DictReader(open((csv_file),'r'))

	for row in counts_data:
		tax_id = row['tax_id']
		# create a dictionary for each taxid if seen for the first time:
		if tax_id not in CountsDict:
			CountsDict[tax_id] = {"S11C1_A_1800":0, "S11C1_C_1800":0, "S14C1_B_2200":0, "S14C1_C_2200":0, "S15C1_B_200":0, "S15C1_C_200":0, "S16C1_A_600":0, "S16C1_B_600":0, "S17C1_A_1000":0, "S17C1_B_1000":0, "S18C1_A_1400":0, "S18C1_C_1400":0, "S19C1_A_1800":0, "S19C1_C_1800":0, "S20C1_B_2200":0, "S20C1_C_2200":0, "S21C1_B_200":0, "S21C1_C_200":0, "S22C1_A_600":0, "S22C1_C_600":0, "S23C1_B_1000":0, "S23C1_C_1000":0, "S24C1_A_1400":0, "S24C1_B_1400":0, "S26C1_A_1800":0, "S26C1_C_1800":0, "S28C1_B_2200":0, "S28C1_C_2200":0, "S29C1_A_200":0, "S29C1_C_200":0, "S30C1_A_600":0, "S30C1_C_600":0, "S31C1_A_1000":0, "S31C1_C_1000":0, "S32C1_B_1400":0, "S32C1_C_1400":0, "S33C1_A_1800":0, "S33C1_C_1800":0, "S34C1_B_2200":0, "S34C1_C_2200":0, "S35C1_A_200":0, "S35C1_C_200":0, "S06C1_A_600":0, "S06C1_C_600":0, "S07C1_A_1000":0, "S07C1_B_1000":0, "S08C1_B_1400":0, "S08C1_C_1400":0}
		# now, add in the counts for each sample:
		for sample in sample_ids:
			CountsDict[tax_id][sample] += float(row[sample])

# Create a dictionary for each taxid to store counts
# note that S06C1, etc changes to STN ids have been made upstream.
CountsDict = {}
sample_ids = ["S11C1_A_1800","S11C1_C_1800","S14C1_B_2200","S14C1_C_2200","S15C1_B_200","S15C1_C_200","S16C1_A_600","S16C1_B_600","S17C1_A_1000","S17C1_B_1000","S18C1_A_1400","S18C1_C_1400","S19C1_A_1800","S19C1_C_1800","S20C1_B_2200","S20C1_C_2200","S21C1_B_200","S21C1_C_200","S22C1_A_600","S22C1_C_600","S23C1_B_1000","S23C1_C_1000","S24C1_A_1400","S24C1_B_1400","S26C1_A_1800","S26C1_C_1800","S28C1_B_2200","S28C1_C_2200","S29C1_A_200","S29C1_C_200","S30C1_A_600","S30C1_C_600","S31C1_A_1000","S31C1_C_1000","S32C1_B_1400","S32C1_C_1400","S33C1_A_1800","S33C1_C_1800","S34C1_B_2200","S34C1_C_2200","S35C1_A_200","S35C1_C_200","S06C1_A_600","S06C1_C_600","S07C1_A_1000","S07C1_B_1000","S08C1_B_1400","S08C1_C_1400"]

# go through the 24 files and add the counts to the dictionary for each taxid
list_file = open(args.list,'r')
for csv_file in list_file:
	add_file_counts(csv_file.strip())

# after adding the counts, create a total summed field
for tax_id in CountsDict:
	CountsDict[tax_id]["total_counts"] = 0
	for sample in sample_ids:
		CountsDict[tax_id]["total_counts"] += CountsDict[tax_id][sample]

# output:
output_csv = open(args.output, 'w')
# this is the same ordering as sample_ids above, but with tax_id, and S6C1 -> S06C1 fixed, etc
header = ",".join(["tax_id","total_counts","S11C1_A_1800","S11C1_C_1800","S14C1_B_2200","S14C1_C_2200","S15C1_B_200","S15C1_C_200","S16C1_A_600","S16C1_B_600","S17C1_A_1000","S17C1_B_1000","S18C1_A_1400","S18C1_C_1400","S19C1_A_1800","S19C1_C_1800","S20C1_B_2200","S20C1_C_2200","S21C1_B_200","S21C1_C_200","S22C1_A_600","S22C1_C_600","S23C1_B_1000","S23C1_C_1000","S24C1_A_1400","S24C1_B_1400","S26C1_A_1800","S26C1_C_1800","S28C1_B_2200","S28C1_C_2200","S29C1_A_200","S29C1_C_200","S30C1_A_600","S30C1_C_600","S31C1_A_1000","S31C1_C_1000","S32C1_B_1400","S32C1_C_1400","S33C1_A_1800","S33C1_C_1800","S34C1_B_2200","S34C1_C_2200","S35C1_A_200","S35C1_C_200","S06C1_A_600","S06C1_C_600","S07C1_A_1000","S07C1_B_1000","S08C1_B_1400","S08C1_C_1400"])
output_csv.write(header + "\n")

for tax_id in CountsDict.keys():
	outline_elts = [str(tax_id)]
	outline_elts.append(str(CountsDict[tax_id]["total_counts"]))
	for sample in sample_ids:
		outline_elts.append(str(CountsDict[tax_id][sample]))
	output_csv.write(",".join(outline_elts) + "\n")
