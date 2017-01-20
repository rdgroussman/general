#!/usr/bin/env python

# reset_xml.py

# Reset XML files
# simply loads and saves out XML files (resets the XML format from one line to many)

from Bio import Phylo
import sys

def build_output_handle(infile_path):
    handle_elts = infile_path.split(".")
    handle_elts.insert(-1,"r")
    out_xml_path = ".".join(handle_elts)
    return out_xml_path

input_xml_path = sys.argv[1]
out_xml_path = build_output_handle(input_xml_path)

tree = Phylo.read(input_xml_path, 'phyloxml')

Phylo.write(tree, out_xml_path, "phyloxml")
