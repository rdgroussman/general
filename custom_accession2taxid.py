#!/usr/bin/env python

# custom_accession2taxid.py


"""
3/1/2018
Ryan D. Groussman

The purpose is to create a custom 'accession2taxid' file, replicating NCBI's
version of this file linking accession IDs with taxonomy ID's.

Here, we will link the accession of a custom reference database (beginning w/ MarRef_w_virus)
to the taxid values that are included in the defline and output a format like this:

accession       accession.version       taxid   gi
WP_070124587    WP_070124587.1  1656094 1079317924
WP_070124589    WP_070124589.1  1656094 1079317926
WP_070124590    WP_070124590.1  28901   1079317927


According to correspondence with the DIAMOND author Ben Buchfink:
'the standard 'prot.accession2taxid.gz' file has four fields: "accession       accession.version       taxid   gi'
you can totally build your own custom mapping file. Diamond uses the
second and third field in that file, the first and the forth are
ignored.
"""

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("input_txt", help="txt file of FASTA deflines, without the > symbol")
parser.add_argument("-o", "--output", help="Name of output file", type=str)
args = parser.parse_args()

# load a file with the deflines from a reference file, like MarRef_w_virus
deflines = open(args.input_txt, 'r')
outfile = open(args.output, 'w')

# add the header to the output file:
outheader = "accession\taccession.version\ttaxid\tgi\n"
outfile.write(outheader)

# examples:
"""
SAR116_cluster_alpha_proteobacterium_HIMB100_tax909943_locID_2503127802_SAR116_HIMB100_00023760_seqID1
Erythrobacter_litoralis_HTCC2594_tax314225_locID_637857414_ELI_10690_seqID15581
Clostridium_phage_phiCD111_refYP_009208367.1_tax10239
"""

# for each line in the header file
for line in deflines:
	taxid_found = False
	# split by underscore and grab the taxid:
	accession = line.strip()
	line_elts = line.split("_")
	for elt in line_elts:
		if elt.startswith("tax"):
			taxid = elt[3:].strip()
			taxid_found = True
	# let's see if we can get away with an empty first and fourth field:
	if taxid_found == True:
		outline = accession + "\t" + accession + "\t" + taxid + "\t" + "1000000000\n"
		outfile.write(outline)
