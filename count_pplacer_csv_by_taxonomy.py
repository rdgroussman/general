#!/usr/bin/env python

import argparse


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

parser = argparse.ArgumentParser()
parser.add_argument("input_csv", help="CSV-formatted input file")
parser.add_argument("-g", "--high_level_groups", help="Print high-level group results", action="store_true")
parser.add_argument("-t", "--top_taxa", help="Print top taxa per group", action="store_true")
parser.add_argument("-i", "--ingroup_edges", help="Only count edge numbers in this file", type=str)

args = parser.parse_args()
other_counter = 0
outgroup_counter = 0

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

def print_high_level_results():
	"""Prints the high-level taxonomic results to stdout"""

	# Spit out results!
	counts_results = [fix_sample_name(sample_name)]
	# sample_name,
	for group_name in sorted(TaxaCountsDict.keys()):
		norm_count = normalize_counts(summed_sample_name, TaxaCountsDict[group_name])
		norm_count_other = normalize_counts(summed_sample_name, other_counter)
		norm_count_outgroup = normalize_counts(summed_sample_name, outgroup_counter)
		counts_results.append(str(norm_count))
	counts_results.append(str(norm_count_other))
	counts_results.append(str(norm_count_outgroup))
	# line below prints the header
	# print 'sample_name' + "," + ",".join(sorted(TaxaCountsDict.keys())) + ",Other,Outgroups"

	# here's our results...
	print ",".join(counts_results)

def print_top_taxa_per_group(group):
	"""Given a group, will run through all of the tax_id that run under its
	umbrella. Prints out a list of tax_ids belonging to group along with
	the proportion of hits in the file (a quick way to normalize)"""

	# header
	top_taxa_header = "sample_name,group,tax_id,fraction_of_all_hits"

	# for our taxa in TreecolorDict
	for taxa in TreecolorDict[group]:
		if taxa in EdgeCountDict.keys():
			hit_fraction = EdgeCountDict[taxa] / float(line_counter)
			print fix_sample_name(sample_name) + "," + group + "," + str(taxa) + "," + str(hit_fraction)

def add_classification_to_edgedict(classification):
	"""Given an NCBI tax_id classification, will add it to a dictionary if
	it's not there with value 1, otherwise will increment counter"""

	if classification in EdgeCountDict.keys():
		EdgeCountDict[classification] += 1
	elif classification not in EdgeCountDict:
		EdgeCountDict[classification] = 1

def process_ingroup_list(ingroup_file_path):
	""" Process a list of edge numbers for the ingroup and their inclusive
	internal edge numbers, as noted in a .jplace file. If given as
	argument, this script will only count these edges. Example:
	72
	73
	74
	...
	"""

	ingroup_edges_file = open(args.ingroup_edges, 'r')
	ingroup_edges_set = set([])
	for edge in ingroup_edges_file:
		edge = edge.strip()
		ingroup_edges_set.add(edge)
	return ingroup_edges_set

def search_and_count_treecolor_dict(tax_id):
	"""go through each tax list in TreecolorDict
	if it matches a list, increment tally in TaxaCountsDict
	"""
	global other_counter
	global outgroup_counter

	taxid_found = False
	for group_name in TreecolorDict:
		if tax_id in TreecolorDict[group_name]:
			TaxaCountsDict[group_name] += 1
			taxid_found = True
	if taxid_found == False:
		other_counter += 1

def fix_sample_name(sample_name):
	"""Fixes some sample names to standardize format.
	E.g. S8C1 > S08C1, etc """

	if sample_name.startswith("S6C1"):
		sample_name = sample_name.replace("S6C1", "S06C1")
	if sample_name.startswith("S7C1"):
		sample_name = sample_name.replace("S7C1", "S07C1")
	if sample_name.startswith("S8C1"):
		sample_name = sample_name.replace("S8C1", "S08C1")
	return sample_name


# if ingroup list file is provided, load the file:
if args.ingroup_edges != None:
	ingroup_edges_set = process_ingroup_list(args.ingroup_edges)

# parse treecolor info
TreecolorDict = parse_treecolor_info(TCInfo)

# Initialize TaxaCountsDict and Other count and Outgroup count
TaxaCountsDict = initialize_taxa_counts_dict()

# Initialize NormFactorsDict
NormFactorsDict = initialize_norm_counts_dict()

# now open up the guppy csv file
input_csv_path = args.input_csv
in_csv = open(input_csv_path, 'r')

# Counts per edge: this dict stores raw counts for each edge.
EdgeCountDict = {}

summed_sample_name = False
line_counter = 0
in_csv.readline()   # skip the first line
for line in in_csv:
	line_counter += 1
	line_elts = line.split(",")
	classification = line_elts[10]
	edge_num = line_elts[3]
	add_classification_to_edgedict(classification)
	# pull out the origin name (from element [0][1])
	origin = line_elts[0]
	sample_name = origin.split(".")[1]
	if summed_sample_name == False:
		summed_sample_name = sample_name
	elif summed_sample_name != sample_name:
		print "Multiple samples in CSV file!"
		summed_sample_name = "Multiple"
	if args.ingroup_edges != None:
		if edge_num in ingroup_edges_set:
			search_and_count_treecolor_dict(classification)
		elif edge_num not in ingroup_edges_set:
			outgroup_counter += 1
	elif args.ingroup_edges == None:
		search_and_count_treecolor_dict(classification)

if args.high_level_groups == True:
	print_high_level_results()
elif args.top_taxa == True:
	for group in TreecolorDict.keys():
		print_top_taxa_per_group(group)
		# also 'other' and 'outgroup'
	print fix_sample_name(sample_name) + "," + "Other" + "," + "1" + "," + str(other_counter / float(line_counter))
	print fix_sample_name(sample_name) + "," + "Outgroups" + "," + "1" + "," + str(outgroup_counter / float(line_counter))
