#!/usr/bin/env python

# combine_d1_contig_data_FPKMx.py

"""
this script will combine taxonomic and functional annotations
along with FPKMx-normalized kallisto counts for a given set of contigs,
perhaps all grouped by pfam_id

This is our palette:
 PF00504.d1.contig_lengths.csv
 PF00504.d1.rain.csv
 PF00504.d1.kallisto_est_counts.tsv
 PF00504.d1.all_data.csv
 PF00504.d1.normed_counts.csv
 PF00504.d1.nt_contig_ids.txt
 PF00504.d1.lca_tax.tab
 PF00504.d1.contig_ids.txt
 PF00504.d1.hmm_out.tab

We'll also need to 'taxify' using a taxa csv

We want to take these and output one big CSV file with columns like this:

# old style:
contig_id,pfam_id,pfam_name,hmmer_eval,tax_id,tax_name,tax_eval,est_counts1,est_counts2...est_counts[i]

# new output:
contig_id,contig_length,pfam_id,pfam_name,tax_id,tax_name,est_counts1,est_counts2...est_counts[i]

"""

import argparse
import csv


def parse_pfam(ContigData, pfam_data):

	# Pfam headers and example (not included in this file) for a .d1.best_pf.hmm_out.csv file:
	"""
	aa_id,pfam_name,pfam_id,pfam_score,pfam_eval
	S06C1_TRINITY_DN2036211_c0_g3_i1_1,PsaA_PsaB,PF00223.18,146.6,2.6e-40
	S06C1_TRINITY_DN2041570_c1_g1_i4_2,PsaA_PsaB,PF00223.18,162.7,3.5e-45
	S06C1_TRINITY_DN2041570_c1_g2_i1_3,PsaA_PsaB,PF00223.18,177.2,1.4e-49
	S06C1_TRINITY_DN2063517_c0_g1_i3_1,PsaA_PsaB,PF00223.18,107.8,1.3e-28
	"""
	# Extract the wanted data from each line:
	for line in pfam_data:
		line_elts = line.split(",")
		contig_id = line_elts[0]
		pfam_name = line_elts[1]
		pfam_id = line_elts[2]
		pfam_eval = line_elts[4].strip()
		pfam_score = line_elts[3]
		# now populate the ContigData:
		ContigData[contig_id]['pfam_name'] = pfam_name
		ContigData[contig_id]['pfam_id'] = pfam_id
		ContigData[contig_id]['pfam_eval'] = pfam_eval
		ContigData[contig_id]['pfam_score'] = pfam_score

def init_taxa_dict(tax_data):

	TaxDict = {}
	for row in tax_data:
		tax_id = row['tax_id']
		TaxDict[tax_id] = row
	return TaxDict

def init_counts_dict(counts_data):

	CountsDict = {}
	for row in counts_data:
		nt_id = row['target_id']
		CountsDict[nt_id] = row
	return CountsDict

def init_rain_dict(rain_data):

	RainDict = {}
	for row in rain_data:
		nt_id = row['nt_id'].strip('"')
		RainDict[nt_id] = row
	return RainDict

def init_length_dict(length_data):

	LengthDict = {}
	for row in length_data:
		nt_id = row['target_id'].strip('"')
		LengthDict[nt_id] = row
	return LengthDict

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
		# get the taxonomic name, rank, and various phylogenetic information from the taxid:
		if tax_id in TaxDict.keys():
			# we can get these straight-up:
			for field in ['tax_name','tax_rank']:
				ContigData[contig_id][field] = TaxDict[tax_id].get(field, "NA")
			# for these, we'll retrieve their tax_name:
			for field in ['phylum', 'class','order','family','genus']:
				field_id = TaxDict[tax_id].get(field, "NA")
				if len(field_id) > 0:
					ContigData[contig_id][field] = TaxDict[field_id].get("tax_name", "NA")
				else:
					ContigData[contig_id][field] = "NA"
			ContigData[contig_id]['phylo_group'],ContigData[contig_id]['phylo_id'] = check_major_phylo_groups(TaxDict, tax_id)

def taxify_taxids(TaxDict, tax_id):

	if tax_id in TaxDict.keys():
		return TaxDict[tax_id]["tax_name"],TaxDict[tax_id]["rank"]
	else:
		return "NA","NA"

def parse_rain(ContigData, RainDict):

	rain_fields = ["pVal","phase","BH"]
	for contig_id in ContigData:
		nt_id = contig_id[:-2]
		if nt_id in RainDict.keys():
			for field in rain_fields:
				ContigData[contig_id][field] = RainDict[nt_id][field]
			if float(ContigData[contig_id]["pVal"]) < 0.05 and float(ContigData[contig_id]["BH"]) < 0.1:
				ContigData[contig_id]["rain_sig"] = True
			else:
				ContigData[contig_id]["rain_sig"] = False
		elif nt_id not in RainDict.keys():
			ContigData[contig_id]["rain_sig"] = False

def parse_contig_lengths(ContigData, LengthDict):

	for contig_id in ContigData:
		nt_id = contig_id[:-2]
		if nt_id in LengthDict.keys():
			ContigData[contig_id]["contig_length"] = LengthDict[nt_id]["length"]

def check_major_phylo_groups(TaxDict, tax_id):

	# Here, we want to see if this contig has membership in
	# taxonomic groups of particular interest:
	"""
	SQL gives us the answer here:
	SELECT tax_name, tax_id, rank FROM taxonomy WHERE tax_name IN ('Bacillariophyta', 'Dinophyceae', 'Haptophyceae', 'Chlorophyta', 'Pelagophyceae', 'Bacteria', 'Archaea', 'Viruses','Opisthokonta');
	Archaea|2157|superkingdom
	Bacteria|2|superkingdom
	Viruses|10239|superkingdom
	Haptophyceae|2830|below_superkingdom
	Opisthokonta|33154|below_superkingdom
	Bacillariophyta|2836|phylum
	Chlorophyta|3041|phylum
	Dinophyceae|2864|class
	Pelagophyceae|35675|class

	# adding in a few more:
	39119	Dictyochophyceae

	"""

	# these are mutually exclusive groups, so we can check each one then:
	if tax_id in TaxDict.keys():
		if TaxDict[tax_id]['superkingdom'] == '2157':
			return "Archaea","2157"
		elif TaxDict[tax_id]['superkingdom'] == '2':
			return "Bacteria","2"
		elif TaxDict[tax_id]['superkingdom'] == '10239':
			return "Viruses","10239"
		elif TaxDict[tax_id]['below_superkingdom'] == '2830':
			return "Haptophyceae","2830"
		elif TaxDict[tax_id]['phylum'] == '2836':
			return "Bacillariophyta","2836"
		elif TaxDict[tax_id]["class"] == "39119":
			return "Dictyochophyceae","39119"
		elif TaxDict[tax_id]['class'] == '35675':
			return "Pelagophyceae","35675"
		elif TaxDict[tax_id]['phylum'] == '3041':
			return "Chlorophyta","3041"
		elif TaxDict[tax_id]['class'] == '2864':
			return "Dinophyceae","2864"
		elif TaxDict[tax_id]['below_superkingdom'] == '33154':
			return "Opisthokonta","33154"
		elif len(TaxDict[tax_id]['phylum']) > 0:
				return TaxDict[TaxDict[tax_id]['phylum']]['tax_name'],TaxDict[TaxDict[tax_id]['phylum']]['tax_id']
		else:
			return "NA","NA"

def count_stats(CountsDict):

	# Add the min, max, mean of all, and the per-station and per-time stats:
	count_rows = ["S11C1_A_1800","S11C1_C_1800","S14C1_B_2200","S14C1_C_2200","S15C1_B_200","S15C1_C_200","S16C1_A_600","S16C1_B_600","S17C1_A_1000","S17C1_B_1000","S18C1_A_1400","S18C1_C_1400","S19C1_A_1800","S19C1_C_1800","S20C1_B_2200","S20C1_C_2200","S21C1_B_200","S21C1_C_200","S22C1_A_600","S22C1_C_600","S23C1_B_1000","S23C1_C_1000","S24C1_A_1400","S24C1_B_1400","S26C1_A_1800","S26C1_C_1800","S28C1_B_2200","S28C1_C_2200","S29C1_A_200","S29C1_C_200","S30C1_A_600","S30C1_C_600","S31C1_A_1000","S31C1_C_1000","S32C1_B_1400","S32C1_C_1400","S33C1_A_1800","S33C1_C_1800","S34C1_B_2200","S34C1_C_2200","S35C1_A_200","S35C1_C_200","S06C1_A_600","S06C1_C_600","S07C1_A_1000","S07C1_B_1000","S08C1_B_1400","S08C1_C_1400"]

def parse_raw_counts(ContigData, CountsDict):

	# this is hard-coded for now. these are the columns for raw 'est_count' kallisto counts
	# from head -1 /mnt/nfs/ryan/diel1/completed_assemblies/transrate_QCed/clustered/bestframe_id99/kallisto/0200_counts_full/combined/diel1.0200.bf100_id99.kallisto_est_counts.tsv
	# I've changed S6C1 -> S06C1
	norm_cols="target_id,S11C1_A_1800,S11C1_C_1800,S14C1_B_2200,S14C1_C_2200,S15C1_B_200,S15C1_C_200,S16C1_A_600,S16C1_B_600,S17C1_A_1000,S17C1_B_1000,S18C1_A_1400,S18C1_C_1400,S19C1_A_1800,S19C1_C_1800,S20C1_B_2200,S20C1_C_2200,S21C1_B_200,S21C1_C_200,S22C1_A_600,S22C1_C_600,S23C1_B_1000,S23C1_C_1000,S24C1_A_1400,S24C1_B_1400,S26C1_A_1800,S26C1_C_1800,S28C1_B_2200,S28C1_C_2200,S29C1_A_200,S29C1_C_200,S30C1_A_600,S30C1_C_600,S31C1_A_1000,S31C1_C_1000,S32C1_B_1400,S32C1_C_1400,S33C1_A_1800,S33C1_C_1800,S34C1_B_2200,S34C1_C_2200,S35C1_A_200,S35C1_C_200,S06C1_A_600,S06C1_C_600,S07C1_A_1000,S07C1_B_1000,S08C1_B_1400,S08C1_C_1400"
	count_rows = ["S11C1_A_1800","S11C1_C_1800","S14C1_B_2200","S14C1_C_2200","S15C1_B_200","S15C1_C_200","S16C1_A_600","S16C1_B_600","S17C1_A_1000","S17C1_B_1000","S18C1_A_1400","S18C1_C_1400","S19C1_A_1800","S19C1_C_1800","S20C1_B_2200","S20C1_C_2200","S21C1_B_200","S21C1_C_200","S22C1_A_600","S22C1_C_600","S23C1_B_1000","S23C1_C_1000","S24C1_A_1400","S24C1_B_1400","S26C1_A_1800","S26C1_C_1800","S28C1_B_2200","S28C1_C_2200","S29C1_A_200","S29C1_C_200","S30C1_A_600","S30C1_C_600","S31C1_A_1000","S31C1_C_1000","S32C1_B_1400","S32C1_C_1400","S33C1_A_1800","S33C1_C_1800","S34C1_B_2200","S34C1_C_2200","S35C1_A_200","S35C1_C_200","S6C1_A_600","S6C1_C_600","S7C1_A_1000","S7C1_B_1000","S8C1_B_1400","S8C1_C_1400"]


	for contig_id in ContigData:
		nt_id = contig_id[:-2]
		if nt_id in CountsDict.keys():
			for count in count_rows:
				ContigData[contig_id][count] = CountsDict[nt_id][count]
			total_count = sum([float(CountsDict[nt_id][count]) for count in count_rows])
			ContigData[contig_id]['total_count'] = total_count

def fill_gaps(ContigData):

	key_fields = ['contig_length','pfam_score', 'pfam_name', 'pfam_eval', 'pfam_id', 'tax_rank', 'tax_name', 'phylum','class','order','family','genus', 'tax_id', 'go_term', 'total_count', 'rain_sig', 'pVal', 'phase','BH']
	for contig in ContigData:
		for key in key_fields:
			if key not in ContigData[contig]:
				ContigData[contig][key] = "NA"

def fix_station_ids(output_header):

	# the usual, S6C1 > S06C1, etc
	return output_header.replace('S6C1', 'S06C1').replace('S7C1', 'S07C1').replace('S8C1', 'S08C1')

def output_csv(ContigData, out_file_name):

	# for now, we're outputting these fields:
	# contig_id,contig_length,pfam_id,pfam_name,tax_id,tax_name,est_counts1,est_counts2...est_counts[i]
	output_fields = ['contig_length','tax_name', 'tax_id', 'tax_rank','phylum','class','order','family','genus', 'pfam_name', 'pfam_id', 'pfam_score', 'pfam_eval', 'go_term', 'total_count', 'rain_sig', 'pVal', 'phase','BH']
	count_rows = ["S11C1_A_1800","S11C1_C_1800","S14C1_B_2200","S14C1_C_2200","S15C1_B_200","S15C1_C_200","S16C1_A_600","S16C1_B_600","S17C1_A_1000","S17C1_B_1000","S18C1_A_1400","S18C1_C_1400","S19C1_A_1800","S19C1_C_1800","S20C1_B_2200","S20C1_C_2200","S21C1_B_200","S21C1_C_200","S22C1_A_600","S22C1_C_600","S23C1_B_1000","S23C1_C_1000","S24C1_A_1400","S24C1_B_1400","S26C1_A_1800","S26C1_C_1800","S28C1_B_2200","S28C1_C_2200","S29C1_A_200","S29C1_C_200","S30C1_A_600","S30C1_C_600","S31C1_A_1000","S31C1_C_1000","S32C1_B_1400","S32C1_C_1400","S33C1_A_1800","S33C1_C_1800","S34C1_B_2200","S34C1_C_2200","S35C1_A_200","S35C1_C_200","S6C1_A_600","S6C1_C_600","S7C1_A_1000","S7C1_B_1000","S8C1_B_1400","S8C1_C_1400"]
	for count in count_rows:
		output_fields.append(count)
	# initialize our output file:
	output_csv = open(out_file_name, 'w')
	# construct and write out the header:
	output_header = "contig_id" + "," + ",".join(output_fields) + "\n"
	output_header = fix_station_ids(output_header)
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
parser.add_argument("-g", "--go_term", help="Informal name for GO term or functional group", type=str)

args = parser.parse_args()

# capture the GO term if given:
pfam_id = args.pfam_id
if args.go_term == False:
	go_term = "NA"
else:
	go_term = args.go_term

# load in contig_ids
print "Loading contig ids..."
contig_ids = open((pfam_id + ".d1.contig_ids.txt"),'r')

# load in pfam data (hmm_out.tab)
print "Loading pfam data..."
pfam_data = open((pfam_id + ".d1.best_pf.hmm_out.csv"),'r')

# load in diamond data (lca_tax.tab)
print "Loading DIAMOND data..."
dmnd_data = open((pfam_id + ".d1.lca_tax.tab"),'r')

# length data
print "Loading length data..."
length_data = csv.DictReader(open((pfam_id + ".d1.contig_lengths.csv"),'r'))
LengthDict = init_length_dict(length_data)

# load in RAIN data (.d1.rain.csv)
print "Loading RAIN data..."
rain_data = csv.DictReader(open((pfam_id + ".d1.rain.csv"),'r'), fieldnames = ( "nt_id","pVal","phase","peak_shape","period","BH" ))
RainDict = init_rain_dict(rain_data)

# load in raw est. counts and add some summary stats:
print "Loading in raw counts..."
counts_data = csv.DictReader(open((pfam_id + ".d1.kallisto_est_counts.tsv"),'r'), delimiter='\t')
CountsDict = init_counts_dict(counts_data)

# load in taxa.csv file and parse
"Loading and parsing NCBI taxonomy..."
tax_data = csv.DictReader(open(args.taxa_csv))
TaxDict = init_taxa_dict(tax_data)

# Start building a ContigData where we can store data for each contig (with contig_ids as keys)
ContigData = {}

# Create a nested dict for each contig:
print "Creating contig dictionary..."
for contig in contig_ids:
	ContigData[contig.strip()] = {}

# population with PFAM data:
print "Populating with PFAM data..."
parse_pfam(ContigData, pfam_data)

# populate with diamond LCA data and get the taxonomic name:
print "Populating with DIAMOND data..."
parse_dmnd_lca(ContigData, dmnd_data)

# populat with RAIN data:
print "Populating with RAIN data..."
parse_rain(ContigData, RainDict)

# add length data:
print "Populating with contig length data..."
parse_contig_lengths(ContigData, LengthDict)

# combine the counts data with the contig data, including some summary stats:
print "Adding in raw count data..."
parse_raw_counts(ContigData, CountsDict)

# add in the GO term:
for contig in ContigData.keys():
	ContigData[contig]['go_term'] = go_term

# fill in the holes in the data
# e.g., no tax placement given
print "Filling in holes..."
fill_gaps(ContigData)

# now, write out the results to the output file!
print "Preparing to write out results..."
output_csv(ContigData, args.out_file)

print "Finished!"
