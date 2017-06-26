#!/usr/bin/env python


# parse_rapsearch_output.py


"""
Ryan Groussman
Armbrust Lab, April 2017

Purpose: Take as in put a '.m8' file from RAPsearch2 output.

For each 'Query', will keep the best matching 'Subject' hit based on bit-score.
Will keep Query, Subject, log(e-value) and bit-score fields.
Furthermore, will parse 'tax_id' from subject.

Will output these values in CSV format.


"""
import re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-t", "--taxa_csv", help="Input taxa.csv", type=str)
parser.add_argument("-m", "--m8", help="Input .m8 file", type=str)
parser.add_argument("-o", "--out_file", help="Specify the name of the outfile", type=str)
args = parser.parse_args()

# m8_file = open(args.m8, 'r')

kept_fields = ["phylum","class"]

def parse_m8_line(m8_line):
	"""parses one line from a rapsearch m8 file; returns query, subject id, e-value and bit-score as a list"""
	m8_elts = m8_line.split('\t')
	return [m8_elts[0], m8_elts[1].replace(',',''), m8_elts[10], m8_elts[11]]

def get_taxid_from_subject(subject_id):
	tax_pattern = re.compile("(tax[0-9]+)")
	subj_elts = subject_id.split("_")
	for elt in subj_elts:
		if tax_pattern.match(elt):
			tax_id = tax_pattern.findall(elt)[0][3:]
			return tax_id


def parse_taxacsv_header(header):

	column_data = []
	header_elts = header.split(",")
	for i in range(len(header_elts)):
		elt = header_elts[i].strip('"')
		if elt in kept_fields:
			column_data.append((i, elt))
	return column_data

TaxDict = {}
with open(args.taxa_csv, 'r') as taxa_csv:
	header = next(taxa_csv)
	column_data = parse_taxacsv_header(header)
	for line in taxa_csv:
		line_elts = line.split(",")
		tax_id = line_elts[0].strip('""')
		rank = line_elts[2].strip('""')
		tax_name = line_elts[3].strip('""')
		TaxDict[tax_id] = {}
		TaxDict[tax_id]["tax_name"] = tax_name
		TaxDict[tax_id]["tax_id"] = tax_id
		TaxDict[tax_id]["rank"] = rank
		for i in range(len(column_data)):
			TaxDict[tax_id][column_data[i][1].strip('""')] = line_elts[int(column_data[i][0])].strip('""')

m8_dict = {}
with open(args.m8, 'r') as m8:
	for line in m8:
		if line.startswith('#'):
			pass
		else:
			m8_elts = parse_m8_line(line)
			query_id = m8_elts[0]
			subject_id = m8_elts[1]
			evalue = m8_elts[2]
			bitscore = m8_elts[3]
			tax_id = get_taxid_from_subject(subject_id)
			if query_id not in m8_dict.keys():
				m8_dict[query_id] = {}
				m8_dict[query_id]["subject_id"] = subject_id
				m8_dict[query_id]["tax_id"] = tax_id
				m8_dict[query_id]["bitscore"] = float(bitscore)
				m8_dict[query_id]["e-value"] = evalue
			elif query_id in m8_dict.keys():
				if float(bitscore) > m8_dict[query_id]["bitscore"]:
					m8_dict[query_id]["subject_id"] = subject_id
					m8_dict[query_id]["tax_id"] = tax_id
					m8_dict[query_id]["e-value"] = evalue
					m8_dict[query_id]["bitscore"] = float(bitscore)

# print it out:
out_file = open(args.out_file, 'w')
header = "query_id,subject_id,tax_id,bitscore,e-value\n"
out_file.write(header)
for query_id in m8_dict:
	outline = query_id +","+ m8_dict[query_id]["subject_id"] +","+ m8_dict[query_id]["tax_id"] +","+ str(m8_dict[query_id]["bitscore"]) +","+ m8_dict[query_id]["e-value"] + "\n"
	out_file.write(outline)
