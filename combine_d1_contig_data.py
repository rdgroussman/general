#!/usr/bin/env python

# combine_contig_data.py

"""
this script will combine taxonomic and functional annotations
along with normalized kallisto counts for a given set of contigs,
perhaps all grouped by pfam_id

This is our palette:
PF00862.d1.normed_counts.csv
PF00862.d1.nt_contig_ids.txt
PF00862.d1.lca_tax.tab
PF00862.d1.contig_ids.txt
PF00862.d1.hmm_out.tab

We'll also need to 'taxify' using a taxa csv

We want to take these and output one big CSV file with columns like this:

contig_id,pfam_id,pfam_name,hmmer_eval,tax_id,tax_name,tax_eval,normed_counts1,normed_counts2...normed_counts[i]

"""

import argparse

def parse_taxacsv_header(header):
	kept_fields = ["phylum","class"]
	column_data = []
	header_elts = header.split(",")
	for i in range(len(header_elts)):
		elt = header_elts[i].strip('"')
		column_data.append((i, elt))
	return column_data

def init_taxa_dict():


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
	return TaxDict

def parse_pfam(ContigData, pfam_data):

	# Pfam headers and example (not included in this file):

	"""
	#                                                                             --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----
	# target name                      accession  query name           accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
	#              ------------------- ---------- -------------------- ---------- --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------
	S15C1_TRINITY_DN476627_c0_g1_i1_3 -          Sucrose_synth        PF00862.18     0.046   20.3   0.0     0.055   20.0   0.0   1.0   1   0   0   1   1   1   1 len100
	S17C1_TRINITY_DN2353879_c1_g2_i1_1 -          Sucrose_synth        PF00862.18   2.9e-08   40.4   0.0   3.2e-08   40.2   0.0   1.0   1   0   0   1   1   1   1 len171
	"""
	# Extract the wanted data from each line:
	for line in pfam_data:
		line_elts = line.split()
		contig_id = line_elts[0]
		pfam_name = line_elts[2]
		pfam_id = line_elts[3]
		pfam_eval = line_elts[4]
		pfam_score = line_elts[5]
		# now populate the ContigData:
		ContigData[contig_id]['pfam_name'] = pfam_name
		ContigData[contig_id]['pfam_id'] = pfam_id
		ContigData[contig_id]['pfam_eval'] = pfam_eval
		ContigData[contig_id]['pfam_score'] = pfam_score

def parse_dmnd_lca(ContigData, dmnd_data):

	# There's no headers in the LCA file, but it's like this:
	# contig_id	tax_id	eval

	for line in dmnd_data:
		line_elts = line.split()
		contig_id = line_elts[0]
		tax_id = line_elts[1]
		dmnd_eval = line_elts[2]
		ContigData[contig_id]['tax_id'] = tax_id
		ContigData[contig_id]['dmnd_eval'] = dmnd_eval
		# get the taxonomic name and rank from the taxid:
		tax_name,tax_rank = taxify_taxids(TaxDict, tax_id)
		ContigData[contig_id]['tax_name'] = tax_name
		ContigData[contig_id]['tax_rank'] = tax_rank

def taxify_taxids(TaxDict, tax_id):

	if tax_id in TaxDict.keys():
		return TaxDict[tax_id]["tax_name"],TaxDict[tax_id]["rank"]
	elif tax_id not in TaxDict.keys():
		return "NA","NA"

def parse_normed_counts(ContigData, counts_data):

	# this is hard-coded for now. these are the columns for normalied kallisto counts
	# from head -1 /mnt/nfs/ryan/diel1/completed_assemblies/transrate_QCed/clustered/bestframe_id99/kallisto/0200_counts_full/combined/diel1.0200.bf100_id99.normed_counts.csv
	norm_cols="target_id,S11C1_A_1800,S11C1_C_1800,S14C1_B_2200,S14C1_C_2200,S15C1_B_200,S15C1_C_200,S16C1_A_600,S16C1_B_600,S17C1_A_1000,S17C1_B_1000,S18C1_A_1400,S18C1_C_1400,S19C1_A_1800,S19C1_C_1800,S20C1_B_2200,S20C1_C_2200,S21C1_B_200,S21C1_C_200,S22C1_A_600,S22C1_C_600,S23C1_B_1000,S23C1_C_1000,S24C1_A_1400,S24C1_B_1400,S26C1_A_1800,S26C1_C_1800,S28C1_B_2200,S28C1_C_2200,S29C1_A_200,S29C1_C_200,S30C1_A_600,S30C1_C_600,S31C1_A_1000,S31C1_C_1000,S32C1_B_1400,S32C1_C_1400,S33C1_A_1800,S33C1_C_1800,S34C1_B_2200,S34C1_C_2200,S35C1_A_200,S35C1_C_200,S6C1_A_600,S6C1_C_600,S7C1_A_1000,S7C1_B_1000,S8C1_B_1400,S8C1_C_1400"

def fill_gaps(ContigData):

	key_fields = ['pfam_score', 'pfam_name', 'pfam_eval', 'pfam_id', 'tax_rank', 'tax_name', 'dmnd_eval', 'tax_id']
	for contig in ContigData:
		for key in key_fields:
			if key not in ContigData[contig]:
				ContigData[contig][key] = "NA"

def output_csv(ContigData, out_file_name):

	# for now, we're outputting these fields:
	output_fields = ['tax_name', 'tax_id', 'tax_rank', 'dmnd_eval', 'pfam_name', 'pfam_id', 'pfam_score', 'pfam_eval']
	# initialize our output file:
	output_csv = open(out_file_name, 'w')
	# construct and write out the header:
	output_header = "contig_id,tax_name,tax_id,tax_rank,dmnd_eval,pfam_name,pfam_id,pfam_score,pfam_eval\n"
	output_csv.write(output_header)
	# now we'll go through each contig and output the values
	for contig in ContigData:
		output_elts = [contig]
		for field in output_fields:
			output_elts.append(str(ContigData[contig][field]))
		output_line = ",".join(output_elts) + "\n"
		output_csv.write(output_line)

parser = argparse.ArgumentParser()
parser.add_argument("-t", "--taxa_csv", help="Input taxa.csv", type=str)
parser.add_argument("pfam_id", help="A PFAM id like PF00862", type=str)
parser.add_argument("-o", "--out_file", help="Specify the name of the outfile", type=str)
args = parser.parse_args()

pfam_id = args.pfam_id

# load in contig_ids
contig_ids = open((pfam_id + ".d1.contig_ids.txt"),'r')

# load in pfam data (hmm_out.tab)
pfam_data = open((pfam_id + ".d1.hmm_out.tab"),'r')

# load in diamond data (lca_tax.tab)
dmnd_data = open((pfam_id + ".d1.lca_tax.tab"),'r')

# load in the counts file
counts_data = open((pfam_id + ".d1.normed_counts.csv"),'r')

# load in taxa.csv file and parse
TaxDict = init_taxa_dict()
# Start building a ContigData where we can store data for each contig (with contig_ids as keys)
ContigData = {}
# Create a nested dict for each contig:
for contig in contig_ids:
	ContigData[contig.strip()] = {}

# population with PFAM data:
parse_pfam(ContigData, pfam_data)

# populate with diamond LCA data and get the taxonomic name:
parse_dmnd_lca(ContigData, dmnd_data)

# fill in the holes in the data
# e.g., no tax placement given
fill_gaps(ContigData)

# now, write out the results to the output file!
output_csv(ContigData, args.out_file)

# print ContigData
# for key in ContigData.keys():
# 	print ContigData[key].keys()
# get the nt versions of contig_ids - we'll need these for the count data
