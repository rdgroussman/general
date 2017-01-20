#!/usr/bin/env python

# Ryan Groussman
# Armbrust Lab 2016

"""
combine_guppy_csv.py

Takes two files: a CSV file with norm factors for each sample-run, and a list of 'guppy to_csv'-made CSV files
for each sample.

NORM_FACTOR here:
/Users/rgroussman/data/SCOPE/diel1/diel1_standard_counts_BD.csv

List file should look like this (192 total for full run):
LOV8.S14C1_B_2200.H5C5H.TaxID.csv
LOV8.S14C1_B_2200.HF5MW.TaxID.csv
LOV8.S14C1_B_2200.HF7JC.TaxID.csv
LOV8.S14C1_B_2200.HKVHM.TaxID.csv
LOV8.S14C1_C_2200.H2NVM.TaxID.csv
LOV8.S14C1_C_2200.HF73N.TaxID.csv

Will output a CSV file with the following columns:
edge_num, tax_id, sample-run1, sample-run2, etc (x192)

"""

# Parse arguments


# Load files
# takes the norm counts file + a list of CSV files to process

NORM_COUNTS_PATH="/Users/rgroussman/data/SCOPE/diel1/diel1_NORM_FACTORS.csv"
INPUT_LIST="/Users/rgroussman/data/SCOPE/diel1/single_gene/LOV8.ddp/test_combine/csv_list.txt"
# INPUT_LIST="/Users/rgroussman/data/SCOPE/diel1/single_gene/LOV8.ddp/test_combine/test_csv_list.txt"
outfile_path="test_output.csv"

# on bloom:
# NORM_COUNTS_PATH="/mnt/ryan/diel1/diel1_NORM_FACTORS.csv"
# INPUT_LIST="csv_list.txt"
# outfile_path="combined_norm_factors_output.csv"


def process_norm_counts(NORM_COUNTS_PATH):
    """
    Takes in a CSV file with sample_name,NORM_FACTOR as columns:
        sample_name,NORM_FACTOR
        S11C1_A_1800.H5C5H,11297.32
        S11C1_A_1800.HF5MW,8640.68
        S11C1_A_1800.HF7JC,9243.03
        S11C1_A_1800.HKVHM,9490.61

        Returns a dictionary with sample_name as key, dict as value
        sample_name dict has machine runs as keys, norm_factor as value
    """

    NormCountsDict = {}
    norm_counts_csv = open(NORM_COUNTS_PATH, 'r')

    # skip the first line
    norm_counts_csv.readline()
    for line in norm_counts_csv:
        line_elts = line.split(",")
        full_sample_name = line_elts[0].split(".")
        sample_name = full_sample_name[0]
        machine_run = full_sample_name[1]
        norm_factor = line_elts[1].strip()

        # If this is the first incidence of a sample name, create a dictionary for it within the dict:
        if sample_name not in NormCountsDict:
            NormCountsDict[sample_name] = {}
        # Then add the machine run and norm factor:
        NormCountsDict[sample_name][machine_run] = norm_factor

    return NormCountsDict

def process_csv_list(INPUT_LIST):
    """Takes a list of CSV files that looks like this (note, no header):

        LOV8.S14C1_B_2200.H5C5H.TaxID.csv
        LOV8.S14C1_B_2200.HF5MW.TaxID.csv
        LOV8.S14C1_B_2200.HF7JC.TaxID.csv
        LOV8.S14C1_B_2200.HKVHM.TaxID.csv
        LOV8.S14C1_C_2200.H2NVM.TaxID.csv
        LOV8.S14C1_C_2200.HF73N.TaxID.csv
        LOV8.S14C1_C_2200.HFHN2.TaxID.csv
        LOV8.S14C1_C_2200.HL2CJ.TaxID.csv

    Will create and return a dictionary with the full filename as key,
    and dictionary values for sample_name and machine_run:

    Example:

    """

    CSV_FileDict = {}
    csv_file_list = open(INPUT_LIST, 'r')

    for line in csv_file_list:
        filename, sample_name, machine_run = process_sample_run_handle(line)
        # line = line.strip()
        CSV_FileDict[filename] = {}
        # line_elts = line.split(".")
        # sample_name = line_elts[1]
        # machine_run = line_elts[2]
        CSV_FileDict[filename]["sample_name"] = sample_name
        CSV_FileDict[filename]["machine_run"] = machine_run

    return CSV_FileDict

def process_guppy_csv(guppy_csv):
    """
    Actually opens up each of the provided CSV files made by guppy
    Will count the incidence of each edge number, along with the
    tax_id that it belongs to. Writes these values to ReadCountsDict like so:

        ReadCountsDict[sample_name][machine_run]{'edge_num':edge_num, 'tax_id':tax_id, 'raw_count':raw_count, 'norm_count':norm_count }

    """

    # modify ReadCountsDict to include this sample's sample_name and machine_run
    guppy_csv, sample_name, machine_run = process_sample_run_handle(guppy_csv)
    add_guppy_csv_to_dict(sample_name, machine_run)

    # open up the guppy csv
    csv_file = open(guppy_csv, 'r')
    csv_file.readline() # skip first line

    # Go through each line of the CSV file and collect edge_num, tax_id
    for line in csv_file:
        line_elts = line.split(",")
        edge_num = line_elts[3]
        tax_id = line_elts[10]

        # Add these values to dictionary if they're not already there and start counting:
        if edge_num not in ReadCountsDict[sample_name][machine_run]:
            ReadCountsDict[sample_name][machine_run][edge_num] = {}
            ReadCountsDict[sample_name][machine_run][edge_num]['tax_id'] = tax_id
            ReadCountsDict[sample_name][machine_run][edge_num]['raw_count'] = 1

            # and add it to the set of edge_num:
            edge_num_set.add(edge_num)

        # And if the edge_num is already there, increment the counter
        elif edge_num in ReadCountsDict[sample_name][machine_run]:
            ReadCountsDict[sample_name][machine_run][edge_num]['raw_count'] += 1

    # When we've gone through all the lines, we can multiply the raw_count by the norm_factor
    norm_factor = get_norm_factor(sample_name, machine_run)

    for edge_num in ReadCountsDict[sample_name][machine_run]:
        ReadCountsDict[sample_name][machine_run][edge_num]['norm_count'] = int(ReadCountsDict[sample_name][machine_run][edge_num]['raw_count']) * float(norm_factor)

def add_guppy_csv_to_dict(sample_name, machine_run):
    """
    Create nested dictionaries in ReadCountsDict for the sample_name and machine_run

    """

    if sample_name not in ReadCountsDict:
        ReadCountsDict[sample_name] = {}

    if machine_run not in ReadCountsDict[sample_name]:
        ReadCountsDict[sample_name][machine_run] = {}

def process_sample_run_handle(line):
    """Extracts sample name and run_id from a filename"""
    line = line.strip()
    line_elts = line.split(".")
    sample_name = line_elts[1]
    machine_run = line_elts[2]

    return line, sample_name, machine_run

def get_norm_factor(sample_name, machine_run):
    """
    Retrieve the norm factor from NormCountsDict
    """

    return NormCountsDict[sample_name][machine_run]

def get_edge_data(edge_num, sample_run_headers):
    """
    For a given edge_num, retrieves the normalized counts and tax_id for every
    sample and run. Averages the normalized counts between machine runs.

    """

    avg_norm_count_list = []
    for sample_name in sample_run_headers:
        norm_count_sum = 0
        rep_count = 0
        for machine_run in ReadCountsDict[sample_name]:
            rep_count += 1
            if edge_num in ReadCountsDict[sample_name][machine_run]:
                norm_count_sum += ReadCountsDict[sample_name][machine_run][edge_num]['norm_count']
                tax_id = ReadCountsDict[sample_name][machine_run][edge_num]['tax_id']
        if rep_count >= 1:
            avg_norm_count = norm_count_sum / float(rep_count)
        elif rep_count == 0:
            avg_norm_count = 0
        avg_norm_count_list.append(str(avg_norm_count))
    return avg_norm_count_list, tax_id

def write_out_csv(outfile_path):
    """Writes out the data from ReadCountsDict in a CSV file

    Will output a CSV file with the following columns:
    edge_num, tax_id, avg norm counts sample1, avg norm counts sample2, etc
    """

    output_csv = open(outfile_path, 'w')

    # sample_run_headers = []
    sample_headers = []
    for sample_name in ReadCountsDict.keys():
        sample_headers.append(sample_name)
        # for machine_run in ReadCountsDict[sample_name]:
        #     sample_run = ".".join([sample_name, machine_run])
        #     sample_run_headers.append(sample_run)

    # sort sample headers
    sorted_sample_headers = sorted(sample_headers)

    outfile_header = "edge_num,tax_id," + ",".join(sorted_sample_headers) + "\n"

    output_csv.write(outfile_header)

    tax_id = None
    for edge_num in sorted(edge_num_set):
        avg_norm_count_list, tax_id = get_edge_data(edge_num, sorted_sample_headers)
        outline = str(edge_num) + "," + str(tax_id) + "," + ",".join(avg_norm_count_list) + "\n"
        output_csv.write(outline)




def main():

    global edge_num_set
    edge_num_set = set([])

    global ReadCountsDict
    ReadCountsDict = {}

    global NormCountsDict
    NormCountsDict = process_norm_counts(NORM_COUNTS_PATH)
    CSV_FileDict = process_csv_list(INPUT_LIST)

    for csv_file in CSV_FileDict:
        process_guppy_csv(csv_file)

    # # for testing
    # for sample_name in ReadCountsDict:
    #     for machine_run in ReadCountsDict[sample_name]:
    #         for edge in ReadCountsDict[sample_name][machine_run]:
    #             print
    #             print edge
    #             print ReadCountsDict[sample_name][machine_run][edge]

    write_out_csv(outfile_path)

    # print "Number of edges with at least one pquery:"
    # print len(edge_num_set)
    # print "Number of pqueries placed:"
    # pquery_sum = 0
    # for sample_name in ReadCountsDict:
    #     for machine_run in ReadCountsDict[sample_name]:
    #         for edge_num in ReadCountsDict[sample_name][machine_run]:
    #             pquery_sum += int(ReadCountsDict[sample_name][machine_run][edge_num]['raw_count'])
    # print pquery_sum

if __name__ == "__main__":
    main()
