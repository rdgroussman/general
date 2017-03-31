#! /usr/bin/env python

# This script will color the branches in an XML tree file based on NCBI taxonomy for the CAMERA 360 transcriptome data set
# treecolor_diatom_CAM360.py colors designated diatom leaves by their class: 
# Coscinodiscophyceae: blue
# Mediophyceae: green
# Fragilariophyceae: orange
# Bacillariophyceae: red
# Call the script as follows:
# treecolor_diatom_CAM360.py CoolTree.xml
# Output: CoolTree.dcolor.xml

import sys
import re

XMLName = sys.argv[1]

# load tab delimited file containing list and color information
TCInfoPath = "/share/projects/diatom_est/annotation/scripts/phylo_lists/CAM360/treecolor_diatom.csv" 
TCInfo = open(TCInfoPath, 'r')

TCInfoDict = {} # treecolor information dictionary

for line in TCInfo:
	linesplit = line.split("\t")	# tab delimited file
	ListFile = linesplit[0]	# the group name
	ListFilePath = r'/share/projects/diatom_est/annotation/scripts/phylo_lists/CAM360/' + ListFile.strip()
	GroupName = ListFile.split(".")	# pop off the prefix
	NCGRFile = open(ListFilePath, 'r') # open the list file
	NCGRList = []
	for item in NCGRFile: # create a list of MMETSP# from each list file
		item = item.strip()
		NCGRList.append(item)
	GroupName = GroupName[0] # Name of group
	RedVal = linesplit[1] # Collect RGB color values.  
	GreenVal = linesplit[2]
	BlueVal = linesplit[3]
	ColorString = r"<color><red>" + RedVal.strip() + "</red><green>" + GreenVal.strip() + "</green><blue>" + BlueVal.strip() + "</blue></color>" # build ColorString
	TCInfoDict[GroupName] = (ColorString,NCGRList) # Color values and MMETSP# assigned as values to group names
	

# load infiles
InXML = open(XMLName, 'r')

# load outfile
OutString = r"\1." + "dcolor.xml"
OutFileName = re.sub(r"(.+)\.xml",OutString,XMLName)
OutFile = open(OutFileName,'w')


xmlLineCount = 0	# count the progression of lines in the xml
InsertColorHere = 0 # this will be used to insert the color tag

for line in InXML:
	if xmlLineCount not in {0,1,2,3}: #skip the first four lines; internal file name also uses <name>
		if line.find("<name>") >= 0: # if <name> is found on the line
			GSID = "Null"
			LineElements = line.split("_") # separate the ID by underscore
			for item in LineElements:
				if item.find("MMETSP") >= 0: # Find and collect the grindstone id (MMETSP####)
					GSID = item.strip()
			for Group in TCInfoDict.keys(): # search through the tc info dict for presence of GSID
				GroupList = TCInfoDict[Group]
				if GSID in GroupList[1]: # if the GSID is in a list belonging to one of the groups, return the color value for that group
					ColorString = GroupList[0]
					InsertColorHere = xmlLineCount + 2
		elif line.find("<branch_length>") == -1 and InsertColorHere > xmlLineCount: # if branch length is not indicated we'll need to insert the color string earlier
			InsertColorHere -= 1
				
		if xmlLineCount == InsertColorHere:
			OutFile.write(ColorString)
				
	OutFile.write(line)
	xmlLineCount += 1 # increment the loop counter by one
	
	
# # close it all!
InXML.close()
OutFile.close()
