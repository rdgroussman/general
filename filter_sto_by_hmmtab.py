#!/usr/bin/env python

# filter_sto_by_hmmtab.py
# using hmm_out table, retrieve sequences from sto file above a cutoff threshold

import argparse
from Bio import SearchIO
from Bio import SeqIO
from Bio import AlignIO
# format name: hmmer3-tab

def build_output_handle(infile_path):
	suffix = 'fasta'
	if cutoff != None:
		suffix = cutoff[1] + str(int(cutoff[0]))
	handle_elts = infile_path.split(".")
	handle_elts.insert(-1,suffix)
	handle_elts[-1] = "fasta"
	outfile_path = ".".join(handle_elts)
	return outfile_path

parser = argparse.ArgumentParser()
parser.add_argument("input_tab", help="Specify input table from hmmsearch")
parser.add_argument("input_sto", help="Stockholm-format alignment from hmmsearch")
parser.add_argument("-b", "--bitscore", help="Filter by bitscore", type=float)
parser.add_argument("-e", "--evalue", help="Filter by evalue", type=float)
parser.add_argument("-o", "--outfile", help="Specify output filename", type=str)
args = parser.parse_args()

cutoff = None
if args.evalue != None:
	cutoff = (args.evalue, 'eval')
elif args.bitscore != None:
	cutoff = (args.bitscore, 'bitscore')

if args.outfile != None:
	outfile_path = args.outfile
elif args.outfile == None:
	outfile_path = build_output_handle(args.input_sto)
output_sto = open(outfile_path, 'w')

kept_seqs = set([])
hmm_table = SearchIO.read(args.input_tab, 'hmmer3-tab')
for hit in hmm_table:
	if cutoff[1] == 'eval' and hit.evalue <= cutoff[0]:
		kept_seqs.add(hit.id)
	elif cutoff[1] == 'bitscore' and hit.bitscore >= cutoff[0]:
		kept_seqs.add(hit.id)

kept_seqs_list=[]
# load the default stockholm file
for record in SeqIO.parse(args.input_sto, "stockholm"):
	seq_id = record.id.split("/")[0] # This grabs the id before the slash
	if seq_id in kept_seqs:
		kept_seqs_list.append(record)

SeqIO.write(kept_seqs_list, output_sto, "fasta")
output_sto.close()
