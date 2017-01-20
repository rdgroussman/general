#!/usr/bin/env python

import sys


# count_pplacer_csv_by_taxonomy.py


"""
Ryan Groussman
Armbrust Lab, University of Washington 2017

This script will take as input the 'guppy to_csv' CSV output from pplacer jplace
files. Going through the list of pqueries (at the moment designed for the 'keep
at most 1' option) and will count the incidence in each 'classification' field
against a list of NCBI tax_ids for groups of particular taxonomy (e.g. Bacillariophyta)

"""

# load tab delimited file containing list and color information
TCInfoPath = "/Users/rgroussman/Dropbox/Armbrust/bioinfo/scripts/treecolor/MarineRef2_plus_internal/treecolors.csv"
TCInfo = open(TCInfoPath, 'r')

def parse_treecolor_info(TCInfo_CSV):

	TreecolorDict = {} # treecolor information dictionary

	for line in TCInfo:
		line_elts = line.split(",")
		ListFile = line_elts[0]	# the group name
		ListFilePath = r'/Users/rgroussman/Dropbox/Armbrust/bioinfo/scripts/treecolor/MarineRef2_plus_internal/' + ListFile.strip()
		group_name = ListFile.split(".")[0]	# pop off the prefix, this is the clade name
		tax_file = open(ListFilePath, 'r') # open the list file
		tax_list = []
		for tax_id in tax_file:
			tax_id = tax_id.strip()
			tax_list.append(tax_id)
		TreecolorDict[group_name] = tax_list

	return TreecolorDict

def initialize_norm_counts_dict():
	"""from
	diel1_standard_counts_SUMS.txt (via Bryn Durham, 13Jan17)
	Norm factor is the sum of BPDSmall1, BPDSmall2, BPDMedium1, BPDMedium2,
	BPDLarge1 and BPDLarge2 from bowtie2 counts, then averaged, and used to
	divide total 'added' counts along with volume:
	NORM FACTOR = ADDED / AVG / VOLUME
	"""

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

def normalize_counts(sample, raw_sum_count):
	"""
	Normalizes against the combined machine-run sums of counts from
	diel1_standard_counts_SUMS.txt (via Bryn Durham, 13Jan17)
	Norm factor is the sum of BPDSmall1, BPDSmall2, BPDMedium1, BPDMedium2,
	BPDLarge1 and BPDLarge2 from bowtie2 counts, then averaged, and used to
	divide total 'added' counts along with volume:
	NORM FACTOR = ADDED / AVG / VOLUME

	Takes as input the sample_id and raw_count, returns the normalized count.
	"""

	if sample in NormFactorsDict.keys():
		norm_count = raw_sum_count * NormFactorsDict[sample]
	else:
		print "Sample ID not found!"
	return norm_count

def initialize_taxa_counts_dict():
	"""Initialize TaxaCountsDict starting with 0 for desired taxons.
	"""

	# Initialize the taxa counts dictionary
	TaxaCountsDict = {
	'Alveolata' : 0,
	'Amoebozoa' : 0,
	'Bacillariophyta' : 0,
	'Cryptophyta' : 0,
	'Dinophyceae' : 0,
	'Excavata' : 0,
	'Glaucophyta' : 0,
	'Haptophyta' : 0,
	'NonDiatomStramenopile' : 0,
	'NonDinoAlveolata' : 0,
	'Opisthokonta' : 0,
	'Rhizaria' : 0,
	'Rhodophyta' : 0,
	'Viridiplantae' : 0
	}

	return TaxaCountsDict

# parse treecolor info
TreecolorDict = parse_treecolor_info(TCInfo)

# Initialize TaxaCountsDict
TaxaCountsDict = initialize_taxa_counts_dict()

# Initialize NormFactorsDict
NormFactorsDict = initialize_norm_counts_dict()

# now open up the guppy csv file
input_csv_path = sys.argv[1]
in_csv = open(input_csv_path, 'r')

summed_sample_name = False
in_csv.readline()   # skip the first line
for line in in_csv:
	line_elts = line.split(",")
	# pull out the classification (element 10)
	classification = line_elts[10]
	# pull out the origin name (from element [0][1])
	origin = line_elts[0]
	sample_name = origin.split(".")[1]
	if summed_sample_name == False:
		summed_sample_name = sample_name
	elif summed_sample_name != sample_name:
		print "Multiple samples in CSV file!"
		summed_sample_name = "Multiple"
	# go through each tax list in TreecolorDict
	# if it matches a list, increment tally in TaxaCountsDict
	taxid_found = False
	for group_name in TreecolorDict:
		if classification in TreecolorDict[group_name]:
			TaxaCountsDict[group_name] += 1
			taxid_found = True


# Normalize the results through normalize_counts in a new dictionary

# Spit out results!
# results_order = ['Alveolata', 'Amoebozoa', 'Bacillariophyta', 'Cryptophyta', 'Excavata','Glaucophyta','Haptophyta','NonDiatomStramenopile','Opisthokonta','Rhizaria','Rhodophyta','Viridiplantae']

# print results_order
# print sorted(TaxaCountsDict.keys())

counts_results = [sample_name]
# sample_name,
for group_name in sorted(TaxaCountsDict.keys()):
	norm_count = normalize_counts(summed_sample_name, TaxaCountsDict[group_name])
	counts_results.append(str(norm_count))
# line below prints the header
# print 'sample_name' + "," + ",".join(sorted(TaxaCountsDict.keys()))

# here's our results...
print ",".join(counts_results)
