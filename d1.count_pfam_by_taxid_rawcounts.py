#!/usr/bin/env python



# count_pfam_by_taxid.py

"""
This script will take in a file like d1.combined_pfam_dmnd_rawcounts_lengths.csv:
X,aa_id,pfam_name,pfam_id,pfam_score,pfam_eval,tax_id,dmnd_eval,nt_id,
S11C1_A_1800,S11C1_C_1800,S14C1_B_2200,S14C1_C_2200,S15C1_B_200,S15C1_C_200,
S16C1_A_600,S16C1_B_600,S17C1_A_1000,S17C1_B_1000,S18C1_A_1400,S18C1_C_1400,
S19C1_A_1800,S19C1_C_1800,S20C1_B_2200,S20C1_C_2200,S21C1_B_200,S21C1_C_200,
S22C1_A_600,S22C1_C_600,S23C1_B_1000,S23C1_C_1000,S24C1_A_1400,S24C1_B_1400,
S26C1_A_1800,S26C1_C_1800,S28C1_B_2200,S28C1_C_2200,S29C1_A_200,S29C1_C_200,
S30C1_A_600,S30C1_C_600,S31C1_A_1000,S31C1_C_1000,S32C1_B_1400,S32C1_C_1400,
S33C1_A_1800,S33C1_C_1800,S34C1_B_2200,S34C1_C_2200,S35C1_A_200,S35C1_C_200,
S6C1_A_600,S6C1_C_600,S7C1_A_1000,S7C1_B_1000,S8C1_B_1400,S8C1_C_1400,length

# It will create a dictionary for each tax_id, and within that a dictionary for pfam counts
# and within those, nested dictionaries for diel1 sample times

# example:
CountsDict = { 2864: {"Glyco_hydro_3_C":2, "Peptidase_S24":1 } }

"""

def init_taxa_dict(tax_data):

	TaxDict = {}
	for row in tax_data:
		tax_id = row['tax_id']
		TaxDict[tax_id] = row
	return TaxDict


def check_major_phylo_groups(TaxDict, tax_id):

	# Here, we want to see if this contig has membership in
	# taxonomic groups of particular interest:
	"""
	SQL gives us the answer here:
	SELECT tax_name, tax_id, rank FROM taxonomy WHERE tax_name IN ('Bacillariophyta', 'Dinophyceae', 'Haptophyceae', 'Chlorophyta', 'Pelagophyceae', 'Bacteria', 'Archaea', 'Viruses','Opisthokonta');
	Viruses|10239|superkingdom
	Haptophyceae|2830|below_superkingdom
	Bacillariophyta|2836|phylum
	Chlorophyta|3041|phylum
	Dinophyceae|2864|class
	Pelagophyceae|35675|class

	# adding in a few more:
	39119	Dictyochophyceae

	"""

	# these are mutually exclusive groups, so we can check each one then:
	if tax_id in TaxDict.keys():
		if TaxDict[tax_id]['superkingdom'] == '10239':
			return "Viruses","10239"
		elif TaxDict[tax_id]['below_superkingdom'] == '2830':
			return "Haptophyceae","2830"
		elif TaxDict[tax_id]['phylum'] == '2836':
			return "Bacillariophyta","2836"
		elif TaxDict[tax_id]["class"] == "39119":
			return "Dictyochophyceae","39119"
		elif TaxDict[tax_id]['class'] == '35675':
			return "Pelagophyceae","35675"
		elif TaxDict[tax_id]['phylum'] == '3041':
			return "Chlorophyta","3041"
		elif TaxDict[tax_id]['class'] == '2864':
			return "Dinophyceae","2864"
		else:
			return "NA","NA"
	else:
		return "NA","NA"


import argparse
import csv

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--input", help="Input CSV", type=str)
parser.add_argument("-o", "--output", help="Output file", type=str)
parser.add_argument("-t", "--taxa_csv", help="Input taxa.csv", type=str)

args = parser.parse_args()

# First, a blank dictionary and a set of all observed pfams:
CountsDict = {}
pfam_set = set([])

# step 0: load in taxa.csv file and parse
print "Loading and parsing NCBI taxonomy..."
tax_data = csv.DictReader(open(args.taxa_csv))
TaxDict = init_taxa_dict(tax_data)

# now we'll open up this big old csv:
print "Opening combined CSV file..."
input_csv = csv.DictReader(open((args.input),'r'))

line_counter = 0
say_when = [1000,10000,100000,200000,500000,1000000,1500000,2000000,2500000,3000000,3500000,4000000,4500000,5000000,6000000,7000000]
sample_ids = ["S11C1_A_1800","S11C1_C_1800","S14C1_B_2200","S14C1_C_2200","S15C1_B_200","S15C1_C_200","S16C1_A_600","S16C1_B_600","S17C1_A_1000","S17C1_B_1000","S18C1_A_1400","S18C1_C_1400","S19C1_A_1800","S19C1_C_1800","S20C1_B_2200","S20C1_C_2200","S21C1_B_200","S21C1_C_200","S22C1_A_600","S22C1_C_600","S23C1_B_1000","S23C1_C_1000","S24C1_A_1400","S24C1_B_1400","S26C1_A_1800","S26C1_C_1800","S28C1_B_2200","S28C1_C_2200","S29C1_A_200","S29C1_C_200","S30C1_A_600","S30C1_C_600","S31C1_A_1000","S31C1_C_1000","S32C1_B_1400","S32C1_C_1400","S33C1_A_1800","S33C1_C_1800","S34C1_B_2200","S34C1_C_2200","S35C1_A_200","S35C1_C_200","S6C1_A_600","S6C1_C_600","S7C1_A_1000","S7C1_B_1000","S8C1_B_1400","S8C1_C_1400"]

# and iterate through every line:
for row in input_csv:
	# this will let us know how the script is progressing:
	line_counter += 1
	if line_counter in say_when:
		print "Processed ", line_counter
	tax_id = row["tax_id"]
	# check if its in a curated major group; we'll only analyze/keep if so:
	phylo_group, phylo_id = check_major_phylo_groups(TaxDict, tax_id)
	# and only do more if it's not nothing:
	if phylo_group != "NA":
		pfam_id = row["pfam_id"]
		length = row["length"]
		if tax_id not in CountsDict.keys():
			CountsDict[tax_id] = {}
		CountsDict[tax_id]['phylo_group'] = phylo_group
		if pfam_id not in CountsDict[tax_id]:
			CountsDict[tax_id][pfam_id] = {"num_contigs":0, "S11C1_A_1800":0, "S11C1_C_1800":0, "S14C1_B_2200":0, "S14C1_C_2200":0, "S15C1_B_200":0, "S15C1_C_200":0, "S16C1_A_600":0, "S16C1_B_600":0, "S17C1_A_1000":0, "S17C1_B_1000":0, "S18C1_A_1400":0, "S18C1_C_1400":0, "S19C1_A_1800":0, "S19C1_C_1800":0, "S20C1_B_2200":0, "S20C1_C_2200":0, "S21C1_B_200":0, "S21C1_C_200":0, "S22C1_A_600":0, "S22C1_C_600":0, "S23C1_B_1000":0, "S23C1_C_1000":0, "S24C1_A_1400":0, "S24C1_B_1400":0, "S26C1_A_1800":0, "S26C1_C_1800":0, "S28C1_B_2200":0, "S28C1_C_2200":0, "S29C1_A_200":0, "S29C1_C_200":0, "S30C1_A_600":0, "S30C1_C_600":0, "S31C1_A_1000":0, "S31C1_C_1000":0, "S32C1_B_1400":0, "S32C1_C_1400":0, "S33C1_A_1800":0, "S33C1_C_1800":0, "S34C1_B_2200":0, "S34C1_C_2200":0, "S35C1_A_200":0, "S35C1_C_200":0, "S6C1_A_600":0, "S6C1_C_600":0, "S7C1_A_1000":0, "S7C1_B_1000":0, "S8C1_B_1400":0, "S8C1_C_1400":0}
		CountsDict[tax_id][pfam_id]["num_contigs"] += 1
		for sample in sample_ids:
			CountsDict[tax_id][pfam_id][sample] += float(row[sample]) / (float(length)/1000)


# now we'll open up our output file and build the header:
output_csv = open(args.output, 'w')
header = ",".join(["tax_id","pfam_id","phylo_group", "num_contigs","S11C1_A_1800","S11C1_C_1800","S14C1_B_2200","S14C1_C_2200","S15C1_B_200","S15C1_C_200","S16C1_A_600","S16C1_B_600","S17C1_A_1000","S17C1_B_1000","S18C1_A_1400","S18C1_C_1400","S19C1_A_1800","S19C1_C_1800","S20C1_B_2200","S20C1_C_2200","S21C1_B_200","S21C1_C_200","S22C1_A_600","S22C1_C_600","S23C1_B_1000","S23C1_C_1000","S24C1_A_1400","S24C1_B_1400","S26C1_A_1800","S26C1_C_1800","S28C1_B_2200","S28C1_C_2200","S29C1_A_200","S29C1_C_200","S30C1_A_600","S30C1_C_600","S31C1_A_1000","S31C1_C_1000","S32C1_B_1400","S32C1_C_1400","S33C1_A_1800","S33C1_C_1800","S34C1_B_2200","S34C1_C_2200","S35C1_A_200","S35C1_C_200","S06C1_A_600","S06C1_C_600","S07C1_A_1000","S07C1_B_1000","S08C1_B_1400","S08C1_C_1400"])
output_csv.write(header + "\n")

# and go through the big dictionary and write out the counts:
for tax_id in CountsDict.keys():
	for pfam_id in CountsDict[tax_id].keys():
		if pfam_id == "phylo_group":
			pass
		else:
			outline_elts = [str(tax_id), str(pfam_id), CountsDict[tax_id]['phylo_group'], CountsDict[tax_id][pfam_id]["num_contigs"]]
			for sample in sample_ids:
				outline_elts.append(CountsDict[tax_id][pfam_id][sample])
			outline = ",".join(str(elt) for elt in outline_elts) + "\n"
			output_csv.write(outline)

output_csv.close()
