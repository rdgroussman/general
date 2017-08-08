#!/usr/bin/env python

"""
# logfattree.py

Takes as input a 'fat' xml tree from pplacer/guppy
Modifies the branch widths to log space (base 2 default)
Also multiplies by a linear correction factor to beef up fattiness
"""

import argparse
import math
import xml.etree.ElementTree as ET
from Bio import Phylo

parser = argparse.ArgumentParser()
parser.add_argument("input_xml", help="XML-formatted input tree")
parser.add_argument("-b", "--log_base", help="Log base", type=int, default=2)
parser.add_argument("-l", "--linear_correct", help="Linear correction factor", type=float, default=2.0)
parser.add_argument("-o", "--out_file", help="Output file name", type=str)


args = parser.parse_args()

input_xml_path = args.input_xml
log_base = args.log_base
linear_correction = args.linear_correct

def build_output_handle(infile_path):
	handle_elts = infile_path.split(".")
	handle_elts.insert(-1,"log2")
	out_xml_path = ".".join(handle_elts)
	return out_xml_path

if args.out_file != None:
	out_xml_path = args.out_file
else:
	out_xml_path = build_output_handle(input_xml_path)


in_xml = ET.parse(input_xml_path)
root = in_xml.getroot()

for parent in root.getiterator():
	if parent.tag == '{http://www.phyloxml.org}clade' and parent.find('{http://www.phyloxml.org}width') != None:
		raw_width = parent.find('{http://www.phyloxml.org}width').text
		log_width = math.log(float(raw_width),log_base) * linear_correction
		parent.find('{http://www.phyloxml.org}width').text = str(log_width)


# Write back to a file
in_xml.write(out_xml_path, xml_declaration=True)
#################

# Need to 'reset' the xml for strange formatting reasons:
tree = Phylo.read(out_xml_path, 'phyloxml')
Phylo.write(tree, out_xml_path, "phyloxml")
