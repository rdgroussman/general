#!/usr/bin/env python

# sort_dmnd_lca_by_tax.py

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

def build_output_handle(group, taxified_csv_path):
	handle_elts = taxified_csv_path.split(".")
	handle_elts.insert(-1,group)
	outfile_path = ".".join(handle_elts)
	return outfile_path

def scrape_taxified_csv(group):

	# open the taxified_csv
	taxified_csv = open(args.taxified_csv,'r')

	# open the output file:
	out_file_name = build_output_handle(group, args.taxified_csv)
	out_file = open(out_file_name, 'w')

	# now go through each line in the taxified_csv and output if the tax_id is within this certain group:
	for line in taxified_csv:
		line_elts = line.split(",")
		tax_id = line_elts[0]
		if tax_id in TreecolorDict[group]:
			out_file.write(line)

def scrape_others():

	taxified_csv = open(args.taxified_csv,'r')
	out_file_name = build_output_handle("Others", args.taxified_csv)
	out_file = open(out_file_name, 'w')

	# this is our set of all the taxids for our given treecolors groups.
	all_group_taxids = set([])
	for group in TreecolorDict:
		all_group_taxids = all_group_taxids.union(TreecolorDict[group])

	# now we'll output to the others file if they're NOT in this set:
	for line in taxified_csv:
		line_elts = line.split(",")
		tax_id = line_elts[0]
		if tax_id not in all_group_taxids:
			out_file.write(line)



import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-t", "--taxified_csv", help="Input taxified.csv", type=str)
parser.add_argument("-c", "--treecolor_csv", help="Input treecolors csv", type=str)
# parser.add_argument("-o", "--out_file", help="Specify the name of the outfile", type=str)
args = parser.parse_args()


# load tab delimited file containing list and color information
TCInfoPath = args.treecolor_csv

# parse treecolor info
TCInfo = open(TCInfoPath, 'r')
TreecolorDict = parse_treecolor_info(TCInfo)
# test the dictionary for any overlapping sets of tax_ids:
test_treecolor_dict(TreecolorDict)

# for for each group in Treecolors, go through the input taxified.csv file, and scrape out the tax_ids that fall under this group:
for group in TreecolorDict:
	scrape_taxified_csv(group)
# and we'll go through once more and create an output file for tax_ids with no placement in our given groups ('other')
scrape_others()
