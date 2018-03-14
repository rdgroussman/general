#!/usr/bin/env python


# uniq_id_and_mapping.py

"""
Input: FASTA file with 'friendly names' format including NCBI tax_ids. Ex:
	SAR116_cluster_alpha_proteobacterium_HIMB100_tax909943_locID_2503126789_SAR116_HIMB100_00013630_seqID1000
	Scrippsiella_trochoidea_CCMP3099_tax71861_locID_CAMPEP_0115682454_seqID996751
	Megavirus_chiliensis_refYP_004894666.1_tax10239

Process:
1. Parse the defline for each sequence
	a. It will get a unique identifier (uid#) based on an incremental count_tracker
	b. taxid will be parsed from the defline (taxid#)

2. Outputs:
	a. New fasta file with uid# in defline: >uid#### defline
	b. Mapping file 1: uniq_id to tax_id (tab separated, for DIAMOND input)
	c. Mapping file 2: uniq_id to full defline

"""

def parse_header(seq_record_id):

	# increment the count_tracker:
	global count_tracker
	count_tracker += 1
	uid = "uid" + str(count_tracker)
	# split by underscore and grab the taxid:
	accession = seq_record_id.strip()
	line_elts = seq_record_id.split("_")
	for elt in line_elts:
		if elt.startswith("tax"):
			taxid = elt[3:].strip()
			taxid_found = True

	return uid, taxid, accession

def init_uid2taxid():

	# write out the header:
	outheader = "accession\taccession.version\ttaxid\tgi\n"
	uid2taxid.write(outheader)

def init_uid2defline():

	# write out the header:
	outheader = "uniq_id,full_defline\n"
	uid2defline.write(outheader)

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("input_fasta", help="FASTA-formatted input file with friendly headers")
parser.add_argument("-o", "--output_fasta", help="Name of output FASTA with uniq_ids", type=str)
parser.add_argument("-t", "--uid2taxid", help="Mapping file of uniq_id to tax_id, tab separated", type=str)
parser.add_argument("-c", "--uid2defline", help="Mapping file of uniq_id to original defline, comma separated", type=str)
args = parser.parse_args()

# open files and add some headers:
input_fasta = open(args.input_fasta, 'r')
output_fasta = open(args.output_fasta, 'w')
uid2taxid = open(args.uid2taxid, 'w')
init_uid2taxid()
uid2defline = open(args.uid2defline, 'w')
init_uid2defline()

# count_tracker for counting sequences and assigning uids:
# global count_tracker
count_tracker = 0

# iterate through the input sequence file:
for seq_record in SeqIO.parse(input_fasta, "fasta"):
	# get uid, taxid, and defline from seq_record.id
	uid, taxid, defline = parse_header(seq_record.id)
	# write out sequence with new defline (">uid# defline")
	out_record = SeqRecord(seq_record.seq, id= uid, description=defline)
	SeqIO.write(out_record, output_fasta, "fasta")
	# write out uid and taxid to uid2taxid file:
	uid2taxid_line = "NA" + "\t" + uid + "\t" + taxid + "\t" + "NA\n"
	uid2taxid.write(uid2taxid_line)
	# write out uid and defline to uid2defline file:
	uid2defline_line = uid + "," + defline + "\n"
	uid2defline.write(uid2defline_line)
