#! /usr/bin/env python

import sys
import re

XMLName = sys.argv[1]

DiatomListFile = r'/share/projects/diatom_est/annotation/scripts/phylo_lists/diatom.list'
# Stramenopile list excludes diatoms
NonDiatomStramenopileListFile = r'/share/projects/diatom_est/annotation/scripts/phylo_lists/nondiatom_stramenopile.list'
GreenListFile = r'/share/projects/diatom_est/annotation/scripts/phylo_lists/viridiplantae.list'
RedListFile = r'/share/projects/diatom_est/annotation/scripts/phylo_lists/rhodophyta.list'
AlveolataListFile = r'/share/projects/diatom_est/annotation/scripts/phylo_lists/alveolata.list'
# note that the Chromalveolate excludes stramenopiles
ChromalveolataListFile = r'/share/projects/diatom_est/annotation/scripts/phylo_lists/nonstram_chromalveolata.list'
# expand with stramenopiles, chlorlophytes, dinos, etc.

# load infiles
InXML = open(XMLName, 'r')

# load list files
DiatomList = open(DiatomListFile, 'r')
NonDiatomStramenopileList = open(NonDiatomStramenopileListFile, 'r')
GreenList = open(GreenListFile,'r')
RedList = open(RedListFile,'r')
AlveolataList = open(AlveolataListFile,'r')
ChromalveolataList = open(ChromalveolataListFile,'r')

# declare color settings for each group

# for black background:
# DiatomColorString = r"<color><red>255</red><green>102</green><blue>0</blue></color>"
# NonDiatomStramenopileColorString = r"<color><red>255</red><green>255</green><blue>0</blue></color>"
# GreenColorString = r"<color><red>0</red><green>255</green><blue>0</blue></color>"
# RedColorString = r"<color><red>255</red><green>0</green><blue>0</blue></color>"
# AlveolataColorString = r"<color><red>0</red><green>51</green><blue>255</blue></color>"
# ChromalveolataColorString = r"<color><red>0</red><green>0</green><blue>255</blue></color>"

# for a white background:
DiatomColorString = r"<color><red>255</red><green>102</green><blue>0</blue></color>"
NonDiatomStramenopileColorString = r"<color><red>255</red><green>0</green><blue>153</blue></color>"
GreenColorString = r"<color><red>0</red><green>153</green><blue>0</blue></color>"
RedColorString = r"<color><red>255</red><green>0</green><blue>0</blue></color>"
AlveolataColorString = r"<color><red>0</red><green>0</green><blue>255</blue></color>"
ChromalveolataColorString = r"<color><red>0</red><green>0</green><blue>255</blue></color>"


# expand with stramenopiles, chlorlophytes, dinos, etc.

# load outfile
OutString = r"\1." + "color.xml"
OutFileName = re.sub(r"(.+)\.xml",OutString,XMLName)
OutFile = open(OutFileName,'w')

# create dictionary with genera names as key values
DiatomDict = {}
for Entry in DiatomList:
	CleanedEntry = Entry.strip() # remove whitespace
	if CleanedEntry not in DiatomDict: # if it's not already there...
		DiatomDict[CleanedEntry] = 1 # add the entry to the dictionary w arbitrary value
        
NonDiatomStramenopileDict = {}
for Entry in NonDiatomStramenopileList:
	CleanedEntry = Entry.strip() # remove whitespace
	if CleanedEntry not in NonDiatomStramenopileDict: # if it's not already there...
		NonDiatomStramenopileDict[CleanedEntry] = 1 # add the entry to the dictionary w arbitrary value

GreenDict = {}
for Entry in GreenList:
	CleanedEntry = Entry.strip() # remove whitespace
	if CleanedEntry not in GreenDict: # if it's not already there...
		GreenDict[CleanedEntry] = 1 # add the entry to the dictionary w arbitrary value
       
RedDict = {}
for Entry in RedList:
	CleanedEntry = Entry.strip() # remove whitespace
	if CleanedEntry not in RedDict: # if it's not already there...
		RedDict[CleanedEntry] = 1 # add the entry to the dictionary w arbitrary value

AlveolataDict = {}
for Entry in AlveolataList:
	CleanedEntry = Entry.strip() # remove whitespace
	if CleanedEntry not in AlveolataDict: # if it's not already there...
		AlveolataDict[CleanedEntry] = 1 # add the entry to the dictionary w arbitrary value
		
ChromalveolataDict = {}
for Entry in ChromalveolataList:
	CleanedEntry = Entry.strip() # remove whitespace
	if CleanedEntry not in ChromalveolataDict: # if it's not already there...
		ChromalveolataDict[CleanedEntry] = 1 # add the entry to the dictionary w arbitrary value

        
# repeat this for other lists



xmlLineCount = 0	# count the progression of lines in the xml
InsertColorHere = 0 # this will be used to insert the color tag

for line in InXML:
	if xmlLineCount not in {0,1,2,3}: #skip the first four lines; internal file name also uses <name>
		if line.find("<name>") >= 0: # if <name> is found on the line
			LineElements = line.split("_") # separate the ID by underscore
			RawGenus = LineElements[0] # Capture the genus
			RawGenus = RawGenus.strip() # remove whitespace
			Genus = RawGenus.replace("<name>","") # remove <name>
			if str(Genus) in DiatomDict:
				InsertColorHere = xmlLineCount + 2
				ColorString = DiatomColorString
			if str(Genus) in NonDiatomStramenopileDict:
				InsertColorHere = xmlLineCount + 2
				ColorString = NonDiatomStramenopileColorString
			if str(Genus) in GreenDict:
				InsertColorHere = xmlLineCount + 2
				ColorString = GreenColorString
			if str(Genus) in RedDict:
				InsertColorHere = xmlLineCount + 2
				ColorString = RedColorString
			if str(Genus) in AlveolataDict:
				InsertColorHere = xmlLineCount + 2
				ColorString = AlveolataColorString
			if str(Genus) in ChromalveolataDict:
				InsertColorHere = xmlLineCount + 2
				ColorString = ChromalveolataColorString
		elif line.find("<branch_length>") == -1 and InsertColorHere > xmlLineCount: # if branch length is not indicated we'll need to insert the color string earlier
			InsertColorHere -= 1
				
		if xmlLineCount == InsertColorHere:
			OutFile.write(ColorString)
				
	OutFile.write(line)
	xmlLineCount += 1 # increment the loop counter by one
	
	
	
# close it all!
#in file
#lists
InXML.close()
DiatomList.close()
NonDiatomStramenopileList.close()
GreenList.close()
RedList.close()
AlveolataList.close()
OutFile.close()
