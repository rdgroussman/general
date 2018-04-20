#!/usr/bin/env python


# count_lca_tab_by_taxonomy.py

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("input_tab", help="LCA tab-formatted input file")
parser.add_argument("-c", "--treecolor_csv", help="Specify path to an alternate treecolor.csv-like file.", type=str)
parser.add_argument("-s", "--sort", help="Output sorted LCA results into new files defined by treecolors file", action="store_true")

args = parser.parse_args()

def parse_treecolor_info(TCInfo):

	TreecolorDict = {} # treecolor information dictionary

	for line in TCInfo:
		line_elts = line.split(",")
		ListFile = line_elts[0]	# the group name
		ListFilePath = '/'.join(args.treecolor_csv.split('/')[0:-1]) + '/' + ListFile.strip()
		group_name = ListFile.split(".")[0]	# pop off the prefix, this is the clade name
		tax_file = open(ListFilePath, 'r') # open the list file
		tax_set = set([])
		tax_list = []
		for tax_id in tax_file:
			tax_id = tax_id.strip()
			# tax_list.append(tax_id)
			tax_set.add(tax_id)
		TreecolorDict[group_name] = tax_set

	return TreecolorDict

def test_treecolor_dict(TreecolorDict):
	"""
	Tests whether or not there are overlapping tax_id between the different
	groups. Prints a warning to stdout if yes; no message if not."""

	# Test that there are no overlapping tax_ids between the different groups in TreecolorDict:
	tax_group_list = TreecolorDict.keys()
	while len(tax_group_list) > 0:
		any_group = tax_group_list.pop()
		for other_group in tax_group_list:
			overlap = TreecolorDict[any_group].intersection(TreecolorDict[other_group])
			if len(overlap) > 0:
				print "!!! WARNING !!!"
				print "Overlap between ", any_group, "and", other_group
				print overlap

def initialize_taxa_counts_dict(TCInfo):
	"""Build a blank dictionary like default_taxa_counts_dict
	but using specified into in alternate treecolor.csv-like file."""

	TaxaCountsDict = {}
	for line in open(TCInfoPath, 'r'):
		base_file = line.split(",")[0]
		group_name = base_file.split(".")[0]
		TaxaCountsDict[group_name] = 0
	return TaxaCountsDict

def print_high_level_results():
	"""Prints the high-level taxonomic results to stdout"""

	global other_counter

	# Spit out results!
	counts_results = []
	group_list = sorted(TaxaCountsDict.keys())
	for group_name in group_list:
		counts_results.append(str(TaxaCountsDict[group_name]))
	counts_results.append(str(other_counter))
	# line below prints the header
	print ",".join(group_list) + ",Other"
	print ",".join(counts_results)

def add_classification_to_edgedict(classification):
	"""Given an NCBI tax_id classification, will add it to a dictionary if
	it's not there with value 1, otherwise will increment counter"""

	if classification in EdgeCountDict.keys():
		EdgeCountDict[classification] += 1
	elif classification not in EdgeCountDict.keys():
		EdgeCountDict[classification] = 1
		# and if the 'sort' flag is True, make a file to sort out matching output lines:
		# if args.sort == True:


def search_and_count_treecolor_dict(tax_id):
	"""go through each tax list in TreecolorDict
	if it matches a list, increment tally in TaxaCountsDict
	"""
	global other_counter
	global other_tax_set

	taxid_found = False
	for group_name in TreecolorDict:
		if tax_id in TreecolorDict[group_name]:
			TaxaCountsDict[group_name] += 1
			taxid_found = True
	if taxid_found == False:
		other_counter += 1

# parse treecolor info
TCInfoPath = args.treecolor_csv
TCInfo = open(TCInfoPath, 'r')
TreecolorDict = parse_treecolor_info(TCInfo)
# test the dictionary for any overlapping sets of tax_ids:
test_treecolor_dict(TreecolorDict)


# Initialize TaxaCountsDict and Other count and Outgroup count
TaxaCountsDict = initialize_taxa_counts_dict(TCInfo)

# and if the 'sort' flag is True, make a file to sort out matching output lines:
print TaxaCountsDict


# now open up the m8 csv file
input_csv_path = args.input_tab
in_csv = open(input_csv_path, 'r')

# Counts per edge: this dict stores raw counts for each edge.
EdgeCountDict = {}

line_counter = 0
other_counter = 0

# in_csv.readline()   # skip the first line
for line in in_csv:
	line_counter += 1
	line_elts = line.split()
	classification = line_elts[1]
	add_classification_to_edgedict(classification)
	search_and_count_treecolor_dict(classification)

print_high_level_results()

in_csv.close()
