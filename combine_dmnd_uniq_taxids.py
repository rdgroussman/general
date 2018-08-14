#!/usr/bin/env python


# combine_dmnd_uniq_taxids.py

"""
Takes a LIST of LCA files

d1.all_0200.bf100.id99.vs_MarRef.uniq_taxids.txt
d1.all_0600.bf100.id99.vs_MarRef.uniq_taxids.txt
d1.all_1000.bf100.id99.vs_MarRef.uniq_taxids.txt
d1.all_1400.bf100.id99.vs_MarRef.uniq_taxids.txt
d1.all_1800.bf100.id99.vs_MarRef.uniq_taxids.txt
d1.all_2200.bf100.id99.vs_MarRef.uniq_taxids.txt

which each have the format (count taxid)
	   2244 1
     44 100127
      2 1001994
     32 100272
    191 1003142
  27799 1003176
    188 1004268


"""
# keep a running sum of all the counts:
global count_total
count_total = 0


def parse_uniq_taxid_txt(filename):
	global count_total

	# open the taxid counts file and iterate through each line:
	taxid_counts = open(filename, 'r')
	for line in taxid_counts:
		line_elts = line.split()
		count = line_elts[0].strip()
		tax_id = line_elts[1].strip()
		# on the first time a taxid is found, create taxid as key with count as integer value
		if tax_id not in TaxDict.keys():
			TaxDict[tax_id] = int(count)
			count_total += int(count)
		# next time it is seen, increment the count by as much
		elif tax_id in TaxDict.keys():
			TaxDict[tax_id] += int(count)
			count_total += int(count)

# create a dictionary of taxids
TaxDict = {}


# pass through each uniq_taxids file
time_points = ["0200", "0600", "1000", "1400", "1800", "2200"]
for time in time_points:
	filename = "d1.all_" + time + ".bf100.id99.vs_MarRef.uniq_taxids.txt"
	parse_uniq_taxid_txt(filename)

# output in a similar format to the input
outfile_handle="d1.all_times.bf100.id99.vs_MarRef.uniq_taxids.txt"
outfile = open(outfile_handle,'w')

for taxa in TaxDict:
	outline = str(TaxDict[taxa]) + " " + taxa + "\n"
	outfile.write(outline)

print "Total counts:", count_total
