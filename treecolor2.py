#!/usr/bin/env python

import sys
import re
# import xml.etree.ElementTree as ET
from Bio import Phylo

def build_output_handle(infile_path):
    handle_elts = infile_path.split(".")
    handle_elts.insert(-1,"col")
    out_xml_path = ".".join(handle_elts)
    return out_xml_path


input_xml_path = sys.argv[1]

# load tab delimited file containing list and color information
TCInfoPath = "/Users/rgroussman/Dropbox/Armbrust/bioinfo/scripts/treecolor/MarineRef2/treecolors.csv"
TCInfo = open(TCInfoPath, 'r')

TreecolorDict = {} # treecolor information dictionary

for line in TCInfo:
	line_elts = line.split(",")	# tab delimited file
	ListFile = line_elts[0]	# the group name
	ListFilePath = r'/Users/rgroussman/Dropbox/Armbrust/bioinfo/scripts/treecolor/MarineRef2/' + ListFile.strip()
	group_name = ListFile.split(".")[0]	# pop off the prefix, this is the clade name
	tax_file = open(ListFilePath, 'r') # open the list file
	tax_list = []
	for tax_id in tax_file: # create a list of tax_id from each list file
		tax_id = tax_id.strip()
		tax_list.append(tax_id.strip())
	RedVal = line_elts[2] # Collect RGB color values.
	GreenVal = line_elts[3]
	BlueVal = line_elts[4]
	ColorString = r"<color><red>" + RedVal.strip() + "</red><green>" + GreenVal.strip() + "</green><blue>" + BlueVal.strip() + "</blue></color>" + "\n" # build ColorString
	TreecolorDict[group_name] = (ColorString,tax_list) # Color values and MMETSP# assigned as values to group names


# load infiles
# in_xml = open(input_xml_path, 'r')

# load out_xml
out_xml_path = build_output_handle(input_xml_path)
# out_xml = open(out_xml_path,'w')


# go through the xml file.
# if a line has color, change the color to the value in TreecolorDict if there is one.

# write it out and save it.

## example block ##

### Bio.Phylo test code ####
# tree = Phylo.read(input_xml_path, 'phyloxml')
# for clade in tree.find_clades(name=True):
#     print clade.name
# Phylo.write(tree, out_xml_path, "phyloxml")
############################

### test code ### elementtree
# input_xml_path = "test.xml"
# in_xml = ET.parse(input_xml_path)
# root = in_xml.getroot()
#
# for taxon in root.iter('name'):
#     print taxon.attrib
# for child in root[0][1][1]:
#     print child.tag
# # Write back to a file
# in_xml.write(out_xml_path, xml_declaration=True)
#################

def change_default_color(True):
    pass



### OLD CODE ###
################
xmlLineCount = 0	# count the progression of lines in the xml
InsertColorHere = 0 # this will be used to insert the color tag
for line in in_xml:
	if xmlLineCount not in {0,1,2,3}: #skip the first four lines; internal file name also uses <name>
		if line.find("<name>") >= 0: # if <name> is found on the line begin the hunt!
            InsertColor = False
			tax_id = "Null"
			line_elts = line.split("_") # separate the ID by underscore
			for item in line_elts:
				if item.find("tax") >= 0: # Find and collect the grindstone id (NCBI tax_id)
					tax_id = item[3:].strip()
			for Group in TreecolorDict.keys(): # search through the tc info dict for presence of tax_id
				GroupList = TreecolorDict[Group]
				if tax_id in GroupList[1]: # if the tax_id is in a list belonging to one of the groups, return the color value for that group
                    InsertColor = True
					ColorString = GroupList[0]
					InsertColorHere = xmlLineCount + 2
		elif line.find("<branch_length>") == -1 and InsertColorHere > xmlLineCount: # if branch length is not indicated we'll need to insert the color string earlier
			InsertColorHere -= 1

		if xmlLineCount == InsertColorHere:
			out_xml.write(ColorString)

	out_xml.write(line)
	xmlLineCount += 1 # increment the loop counter by one
### OLD CODE ###
################

# # close it all!
