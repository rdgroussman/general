#!/usr/bin/env python

"""
Ryan Groussman
Armbrust Lab
September 2016

This script takes as output a log of memory use created by 'record_mem_use.sh'.
It will output csv with date, used_mem.

Example of input:
    Mon Sep 12 18:40:37 UTC 2016
    8.4G
    Mon Sep 12 18:45:37 UTC 2016
    8.4G
    Mon Sep 12 18:50:37 UTC 2016
    161G
    Mon Sep 12 18:55:37 UTC 2016
    38G
"""

# import packages
import argparse

# parse incoming argument
parser = argparse.ArgumentParser()
parser.add_argument("mem_log", help="log of memory use")
args = parser.parse_args()


# open up the mem.log
mem_log = open(args.mem_log, 'r')

outfile_path = args.mem_log + ".csv"
output_csv = open(outfile_path, 'w')

date_string = ""
mem_string = ""

for line in mem_log:
    if line[:1].isalpha():
        date_string = line.strip()
        date_string = date_string.strip()
    elif line[:1].isdigit:
        mem_string = line.strip().strip("G")
        outstring = date_string + "," + mem_string + "\n"
        output_csv.write(outstring)

output_csv.close()
