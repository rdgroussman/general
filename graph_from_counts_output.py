#!/usr/bin/env python


# graph_from_counts_output.py
# Purpose: generate a number of figures using count data from $GENE.counts_results.csv,
# the output of count_pplacer_csv_by_taxonomy.py


import matplotlib.pyplot as plt
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("input_csv", help="CSV-formatted input file from count_pplacer_csv_by_taxonomy.py")
parser.add_argument("-d", "--diel_expression", help="Generate plot of diel expression over 24h window", action="store_true")
parser.add_argument("-p", "--print_table", help="Print a CSV table of values used to construct graphs", action="store_true")
parser.add_argument("-t", "--taxa_shown", help="Comma-separated string of high-level taxa to include in plot", type=str)
parser.add_argument("-n", "--not_shown", help="Comma-separated string of high-level taxa to NOT include in plot", type=str)

args = parser.parse_args()

def parse_input_csv():
	input_csv_file = open(args.input_csv, 'r')
	CountsResultsDict = {}
	line_counter = 0
	for line in input_csv_file:
		if line_counter == 0:
			header_elts = line.split(",")
		elif line_counter > 0:
			line_elts = line.split(",")
			sample_name = line_elts[0]
			CountsResultsDict[sample_name] = {}
			for i in range(len(header_elts)):
				CountsResultsDict[sample_name][header_elts[i].strip()] = line_elts[i].strip()
			CountsResultsDict[sample_name]['time'] = get_time(sample_name)
		line_counter += 1
	return CountsResultsDict

def get_time(sample_name):
	if sample_name.endswith("_600"):
		return '0600'
	elif sample_name.endswith("_1000"):
		return '1000'
	elif sample_name.endswith("_1400"):
		return '1400'
	elif sample_name.endswith("_1800"):
		return '1800'
	elif sample_name.endswith("_2200"):
		return '2200'
	elif sample_name.endswith("_200"):
		return '0200'

def plot_diel_expression():
	"""Collect the counts for each time point (e.g. 0600, 1000, 1200, etc)
	Plot the expression for all (or specified) groups in this window. Show
	standard error."""

	# First we build the dictionaries.
	# Initialize dict for each time point. Dictionaries within dictionaries.
	DielExpDict = {
	'0200' : {},
	'0600' : {},
	'1000' : {},
	'1400' : {},
	'1800' : {},
	'2200' : {}
	}

	diel_hour_order = sorted(DielExpDict.keys())
	# A list of the values (taxon groups) in CountsResultsDict
	GroupsSet = set([])

	# Go through the CountsResultsDict and collect the counts for each group per time
	for sample in CountsResultsDict:
		sample_time = get_time(sample)
		# If the subsequent keys are NOT 'time', 'sample_name', then continue:
		for key in CountsResultsDict[sample]:
			if key != 'time' and key != 'sample_name':
				GroupsSet.add(key)
				if key not in DielExpDict[sample_time]:
					DielExpDict[sample_time][key] = [int(CountsResultsDict[sample][key])]
				elif key in DielExpDict[sample_time]:
					DielExpDict[sample_time][key].append(int(CountsResultsDict[sample][key]))


	# Now we can start plotting!
	PlotDict = {}
	legend = []

	x = diel_hour_order
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.set_xlabel('Hour', fontsize = 16)
	ax.set_xticks([int(hour) for hour in diel_hour_order])
	print diel_hour_order
	ax.set_ylabel('Estimated Transcripts L-1', fontsize = 16)
	plt.title("DGAT1 Expression over a 24h cycle")

	top_group_avg = 0
	plot_cutoff = 0.1 # Must be this pct of top_group_avg to plot this group

	# First we can collect all the axis data:
	for group in GroupsSet:
		y_vals = [DielExpDict[hour][group] for hour in diel_hour_order] # gives a list of lists
		y_array = [np.array(vals) for vals in y_vals] # np.std needs this array to calculate stdev
		y_mean = [np.mean(list) for list in y_array]
		y_err = [np.std(list)/np.sqrt(len(list)) for list in y_array] # np.std gives stdev of a population; div by sqrt(n) to get s.e.m.
		group_avg = np.mean(y_mean)

		if group_avg >= top_group_avg:
			top_group_avg = group_avg
		PlotDict[group] = {}
		PlotDict[group]['y_mean'] = y_mean
		PlotDict[group]['y_err'] = y_err
		PlotDict[group]['group_avg'] = group_avg

	# only plot if the expression level crosses a threshold:


	for group in GroupsSet:
		if PlotDict[group]['group_avg'] >= top_group_avg * plot_cutoff:

			legend.append(group)
			ax.errorbar(x, PlotDict[group]['y_mean'], yerr = PlotDict[group]['y_err'])


			# plt.errorbar(x, y_mean, yerr=PlotDict[group]['y_err'])
			# plt.plot(x, PlotDict[group]['y_mean'])

			# error bars:
			# http://stackoverflow.com/questions/26378005/matplotlib-getting-different-colors-in-data-lines-with-error-bars



	# old
	# plt.legend(legend, loc='best')
	# plt.show()
	plt.tight_layout()
	plt.savefig("test.png")
	plt.clf()



def plot_expression_over_time_series():
	"""Plot the expression over the full time series for each specified
	taxonomic group. Show error ranges between duplicates."""

CountsResultsDict = parse_input_csv()
plot_diel_expression()
