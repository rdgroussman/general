#!/usr/bin/env python

# Ryan Groussman
# Armbrust Lab 2016

"""
mORFeus.py

Takes as input a FASTA file with translated input, e.g:
>1_1
IA*RKHWRQRWPR*LLSWRLPSWIWKRRRKE*STSKALQRPAKSRSSS*LLRPRNTRMRR
LLTWRGSRRRSNLSERQSRN

Will output each ORF as its own sequence. Length-cutoff criteria
and Met-start arguments can be provided.
"""

import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def build_output_handle(infile_path, min_pep_length):
    handle_elts = infile_path.split(".")
    handle_mod = "orfs" + str(args.min_pep_length)
    handle_elts.insert(-1,handle_mod)
    outfile_path = ".".join(handle_elts)
    return outfile_path

def get_machine_run(infile_path):
    """Given a FASTA path with the machine run included in format:
        ex: S11C1_A_1800.H5C5H_combined.6tr.fasta
    Returns a string with the machine run:
        ex: H5C5H
    """

    string_elts = infile_path.split(".")
    return string_elts[1][:5]

def process_seq(seq_record):
    """
    Takes a fasta-format translation as output by transeq Ex:
        >3_5
        SRQYLQIGFDCRMATTSLPSL*EGDQIHRLHP*EGCISHTR
        SYVHTMLPKNHYQRTNLPECIHGWHSLHHSR**LE*YPX

        QLCGLLTPGDILIAVNGKSLINGTIHNPVSMERMFSVLKPLSQPLDAEDGQYYAREVRLRLVLDEGRELLRDQKEREERK

    Returns individual peptides, accepting arguments for min_pep_length
    and Met start (following stops). Also modifies defline.
    """

    in_frame = seq_record.seq
    orfs = in_frame.split("*")

    good_orfs_count = 0
    starts_with_met = 0

    # if the whole reading frame is uninterrupted, write out
    if len(orfs) == 1:
        good_orfs_count = 1
        write_out_pep(orfs[0], seq_record.id, 1)

    # otherwise, check for length
    elif len(orfs) > 1:
        for orf in orfs:
            if len(orf) >= args.min_pep_length:
                good_orfs_count += 1
                write_out_pep(orf, seq_record.id, good_orfs_count)
                # if orf[0] == "M":
                #     starts_with_met += 1

def write_out_pep(orf, id, counter):

    """Modifies the defline to include the integer number for output ORFs from each reading frame
    as well as the machine run that it came from."""

    run_id = ""
    if args.machine_run == True:
        run_id = get_machine_run(args.fasta_file) + "_"
    defline = ">" + run_id + id + "_" + str(counter) + "\n"
    output_file.write(defline)
    output_file.write(str(orf) + "\n")

# parse incoming argument
parser = argparse.ArgumentParser()
parser.add_argument("fasta_file", help="Translated sequences in FASTA format")
parser.add_argument("-l", "--min_pep_length", help="Minimum peptide length (integer) to keep", type=int)
parser.add_argument("-m", "--machine_run", help="Add Illumina machine run ID to output deflines", action="store_true")
args = parser.parse_args()

# open up the fasta
fasta_file = open(args.fasta_file, 'r')
output_fasta_handle = build_output_handle(args.fasta_file, args.min_pep_length)
output_file = open(output_fasta_handle, 'w')

# for each sequence in the fasta file;
for seq_record in SeqIO.parse(fasta_file, "fasta"):
    process_seq(seq_record)

fasta_file.close()
output_file.close()
