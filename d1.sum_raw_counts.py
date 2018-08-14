#!/usr/bin/env python

# d1.sum_raw_counts.py

"""
Purpose: Sum *all* raw counts for diel1

We'll want to run this one separately on the set of six timepoint collections to avoid overcomplexity.
Once this finishes, it will be trival to combine the counts for each of the six time points.
"""

import argparse
import csv

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--diamond", help="Input DIAMOND LCA file", type=str)
parser.add_argument("-k", "--kallisto", help="Input Kallisto raw counts file", type=str)
parser.add_argument("-o", "--output", help="Output CSV path", type=str)
args = parser.parse_args()

# hard-coded values:
# print a message when entries parsed gets to these thresholds:
counter_markers = [1000, 10000, 100000, 250000, 500000, 1000000, 1500000, 2000000, 2500000, 3000000]
sample_ids = ["S11C1_A_1800","S11C1_C_1800","S14C1_B_2200","S14C1_C_2200","S15C1_B_200","S15C1_C_200","S16C1_A_600","S16C1_B_600","S17C1_A_1000","S17C1_B_1000","S18C1_A_1400","S18C1_C_1400","S19C1_A_1800","S19C1_C_1800","S20C1_B_2200","S20C1_C_2200","S21C1_B_200","S21C1_C_200","S22C1_A_600","S22C1_C_600","S23C1_B_1000","S23C1_C_1000","S24C1_A_1400","S24C1_B_1400","S26C1_A_1800","S26C1_C_1800","S28C1_B_2200","S28C1_C_2200","S29C1_A_200","S29C1_C_200","S30C1_A_600","S30C1_C_600","S31C1_A_1000","S31C1_C_1000","S32C1_B_1400","S32C1_C_1400","S33C1_A_1800","S33C1_C_1800","S34C1_B_2200","S34C1_C_2200","S35C1_A_200","S35C1_C_200","S6C1_A_600","S6C1_C_600","S7C1_A_1000","S7C1_B_1000","S8C1_B_1400","S8C1_C_1400"]

# step 1: dictionary linking nt_id to tax_id (diamond lca files)
# this will also need to choose the BEST tax_id if there are competing annotations (for difft frames)

nt_id2tax_idDict = {} # for nt_id to tax_id
CountsDict = {} # stores values for tax_id counts

increment_counter = 0
with open(args.diamond, 'r') as dmnd_csv:
	for line in dmnd_csv:
		increment_counter += 1
		if increment_counter in counter_markers:
			print "Diamond entries processed: ", increment_counter
		line_elts = line.split()
		nt_id = line_elts[0][:-2] # get nt_id from aa_id
		tax_id = str(line_elts[1].strip())
		if tax_id not in CountsDict.keys():
			CountsDict[tax_id] = {"S11C1_A_1800":0, "S11C1_C_1800":0, "S14C1_B_2200":0, "S14C1_C_2200":0, "S15C1_B_200":0, "S15C1_C_200":0, "S16C1_A_600":0, "S16C1_B_600":0, "S17C1_A_1000":0, "S17C1_B_1000":0, "S18C1_A_1400":0, "S18C1_C_1400":0, "S19C1_A_1800":0, "S19C1_C_1800":0, "S20C1_B_2200":0, "S20C1_C_2200":0, "S21C1_B_200":0, "S21C1_C_200":0, "S22C1_A_600":0, "S22C1_C_600":0, "S23C1_B_1000":0, "S23C1_C_1000":0, "S24C1_A_1400":0, "S24C1_B_1400":0, "S26C1_A_1800":0, "S26C1_C_1800":0, "S28C1_B_2200":0, "S28C1_C_2200":0, "S29C1_A_200":0, "S29C1_C_200":0, "S30C1_A_600":0, "S30C1_C_600":0, "S31C1_A_1000":0, "S31C1_C_1000":0, "S32C1_B_1400":0, "S32C1_C_1400":0, "S33C1_A_1800":0, "S33C1_C_1800":0, "S34C1_B_2200":0, "S34C1_C_2200":0, "S35C1_A_200":0, "S35C1_C_200":0, "S6C1_A_600":0, "S6C1_C_600":0, "S7C1_A_1000":0, "S7C1_B_1000":0, "S8C1_B_1400":0, "S8C1_C_1400":0}
		evalue = float(line_elts[2].strip()) # python handles scientific notation correctly here
		if nt_id not in nt_id2tax_idDict.keys():
			nt_id2tax_idDict[nt_id] = (tax_id, evalue)
		elif evalue < nt_id2tax_idDict[nt_id][1]: # if it exists, and current evalue is smaller, overwrite with current values
			nt_id2tax_idDict[nt_id] = (tax_id, evalue)


# step 2: establish a dictionary for each tax_id, with the 48 samples & counts as key-value pairs
# this is done within the loop, above
# however, we'll establish a special container for 'NA' taxids:
CountsDict["NA"] = {"S11C1_A_1800":0, "S11C1_C_1800":0, "S14C1_B_2200":0, "S14C1_C_2200":0, "S15C1_B_200":0, "S15C1_C_200":0, "S16C1_A_600":0, "S16C1_B_600":0, "S17C1_A_1000":0, "S17C1_B_1000":0, "S18C1_A_1400":0, "S18C1_C_1400":0, "S19C1_A_1800":0, "S19C1_C_1800":0, "S20C1_B_2200":0, "S20C1_C_2200":0, "S21C1_B_200":0, "S21C1_C_200":0, "S22C1_A_600":0, "S22C1_C_600":0, "S23C1_B_1000":0, "S23C1_C_1000":0, "S24C1_A_1400":0, "S24C1_B_1400":0, "S26C1_A_1800":0, "S26C1_C_1800":0, "S28C1_B_2200":0, "S28C1_C_2200":0, "S29C1_A_200":0, "S29C1_C_200":0, "S30C1_A_600":0, "S30C1_C_600":0, "S31C1_A_1000":0, "S31C1_C_1000":0, "S32C1_B_1400":0, "S32C1_C_1400":0, "S33C1_A_1800":0, "S33C1_C_1800":0, "S34C1_B_2200":0, "S34C1_C_2200":0, "S35C1_A_200":0, "S35C1_C_200":0, "S6C1_A_600":0, "S6C1_C_600":0, "S7C1_A_1000":0, "S7C1_B_1000":0, "S8C1_B_1400":0, "S8C1_C_1400":0}

# step 3: go through counts for each nt_id, and add them to the appropriate dictionary for the linked tax_id
increment_counter = 0
# counts_data = csv.DictReader(open((args.kallisto),'r'), delimiter='\t', fieldnames = ("nt_id","S11C1_A_1800","S11C1_C_1800","S14C1_B_2200","S14C1_C_2200","S15C1_B_200","S15C1_C_200","S16C1_A_600","S16C1_B_600","S17C1_A_1000","S17C1_B_1000","S18C1_A_1400","S18C1_C_1400","S19C1_A_1800","S19C1_C_1800","S20C1_B_2200","S20C1_C_2200","S21C1_B_200","S21C1_C_200","S22C1_A_600","S22C1_C_600","S23C1_B_1000","S23C1_C_1000","S24C1_A_1400","S24C1_B_1400","S26C1_A_1800","S26C1_C_1800","S28C1_B_2200","S28C1_C_2200","S29C1_A_200","S29C1_C_200","S30C1_A_600","S30C1_C_600","S31C1_A_1000","S31C1_C_1000","S32C1_B_1400","S32C1_C_1400","S33C1_A_1800","S33C1_C_1800","S34C1_B_2200","S34C1_C_2200","S35C1_A_200","S35C1_C_200","S6C1_A_600","S6C1_C_600","S7C1_A_1000","S7C1_B_1000","S8C1_B_1400","S8C1_C_1400"))
counts_data = csv.DictReader(open((args.kallisto),'r'), delimiter='\t')

for row in counts_data:
	increment_counter += 1
	if increment_counter in counter_markers:
		print "Kallisto entries processed: ", increment_counter
	nt_id = row['target_id'].strip()
	tax_id = nt_id2tax_idDict.get(nt_id, ("NA",0))
	tax_id = tax_id[0]
	# print tax_id
	for sample in sample_ids:
		CountsDict[tax_id][sample] += float(row[sample])


# step 4: save out results to csv:
# tax_id, S8C1_B_1400, S8C1_C_1400, ... S35C1_C_200
output_csv = open(args.output, 'w')
# this is the same ordering as sample_ids above, but with tax_id, and S6C1 -> S06C1 fixed, etc
header = ",".join(["tax_id","S11C1_A_1800","S11C1_C_1800","S14C1_B_2200","S14C1_C_2200","S15C1_B_200","S15C1_C_200","S16C1_A_600","S16C1_B_600","S17C1_A_1000","S17C1_B_1000","S18C1_A_1400","S18C1_C_1400","S19C1_A_1800","S19C1_C_1800","S20C1_B_2200","S20C1_C_2200","S21C1_B_200","S21C1_C_200","S22C1_A_600","S22C1_C_600","S23C1_B_1000","S23C1_C_1000","S24C1_A_1400","S24C1_B_1400","S26C1_A_1800","S26C1_C_1800","S28C1_B_2200","S28C1_C_2200","S29C1_A_200","S29C1_C_200","S30C1_A_600","S30C1_C_600","S31C1_A_1000","S31C1_C_1000","S32C1_B_1400","S32C1_C_1400","S33C1_A_1800","S33C1_C_1800","S34C1_B_2200","S34C1_C_2200","S35C1_A_200","S35C1_C_200","S06C1_A_600","S06C1_C_600","S07C1_A_1000","S07C1_B_1000","S08C1_B_1400","S08C1_C_1400"])
output_csv.write(header + "\n")

for tax_id in CountsDict.keys():
	outline_elts = [str(tax_id)]
	for sample in sample_ids:
		outline_elts.append(str(CountsDict[tax_id][sample]))
	outline = ",".join(outline_elts) + "\n"
	output_csv.write(outline)

output_csv.close()
