#!/usr/bin/env python

"""
# logfattree.py

Takes as input a 'fat' xml tree from pplacer/guppy
Modifies the branch widths to log space (base 2 default)



"""

import argparse
import math
import xml.etree.ElementTree as ET
from Bio import Phylo
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("input_xml", help="XML-formatted input tree")
parser.add_argument("-o", "--out_file", help="Output file name", action="store_true")

args = parser.parse_args()
input_csv_path = args.input_csv


if args.out_file != None:
	out_xml_path = args.out_file
else:
	out_xml_path = build_output_handle(input_xml_path)


in_xml = ET.parse(args.input_xml)
root = in_xml.getroot()

for parent in root.getiterator():
	if parent.tag == '{http://www.phyloxml.org}clade' and parent.find('{http://www.phyloxml.org}width') != None:
		raw_width = parent.find('{http://www.phyloxml.org}width').text
		log_width = math.log(raw_width,2)
		parent.find('{http://www.phyloxml.org}width').text = log_width


# Write back to a file
in_xml.write(out_xml_path, xml_declaration=True)
#################

# Need to 'reset' the xml for strange formatting reasons:
tree = Phylo.read(out_xml_path, 'phyloxml')
Phylo.write(tree, out_xml_path, "phyloxml")
