#!/usr/bin/env python

"""
Ryan Groussman
Armbrust Lab, University of Washington
September 2016

This script compares RAPsearch2 results against NR and MMETSP databases.
Keeps the best hit from each.
Will output one line for each query contig in csv format:
query_id,nr_match,nr_evalue,mmetsp_match,mmetsp_evalue

Will only keep matches if log(e-value) < -20, report 'NA' if none
"""

## testing ###
# nr_m8_path = "/Users/rgroussman/data/SCOPE/testing/test_assemblies/AWS_SCOPE_4lane_1/len750/transrate_results750-3k/RAPsearch_results/test/test_vs_NR.m8"
# mmetsp_m8_path = "/Users/rgroussman/data/SCOPE/testing/test_assemblies/AWS_SCOPE_4lane_1/len750/transrate_results750-3k/RAPsearch_results/test/test_vs_MMETSP.m8"

# real #
# nr_m8_path = "/Users/rgroussman/data/SCOPE/testing/test_assemblies/AWS_SCOPE_4lane_1/len750/transrate_results750-3k/RAPsearch_results/good.SCOPE-4ln.750-3k_vs_NR.rap.m8"
# mmetsp_m8_path = "/Users/rgroussman/data/SCOPE/testing/test_assemblies/AWS_SCOPE_4lane_1/len750/transrate_results750-3k/RAPsearch_results/good.SCOPE-4ln.750-3k_vs_MMETSP.rap.m8"

# on bloom:
mmetsp_m8_path = "/mnt/ryan/assemblies/RAPsearch_vs_MMETSP/good.SCOPE-4ln.750-3k_vs_MMETSP.rap.m8"
nr_m8_path = "/mnt/ryan/assemblies/RAPsearch_vs_NR/good.SCOPE-4ln.750-3k_vs_NR.rap.m8"
outfile_path = "/mnt/ryan/assemblies/RAPsearch_vs_NR/rapsearch_results.csv"


def parse_m8(m8_file):
    """parses the rapsearch-output m8 file.
    Creates a dictionary keeping the best hit and e-value for each query, outputs dictionary
    """

    m8_dict = {}
    for line in m8_file:
        if line.startswith('#'):
            pass
        else:
            m8_elts = parse_m8_line(line)
            query_id = m8_elts[0]
            subject_id = m8_elts[1]
            evalue = m8_elts[2]
            if query_id not in m8_dict.keys():
                m8_dict[query_id] = {}
                m8_dict[query_id]["subject_id"] = subject_id
                m8_dict[query_id]["e-value"] = evalue
            elif query_id in m8_dict.keys():
                if float(evalue) < float(m8_dict[query_id]["e-value"]):
                    m8_dict[query_id]["subject_id"] = subject_id
                    m8_dict[query_id]["e-value"] = evalue
    return m8_dict

def parse_m8_line(m8_line):
    """parses one line from a rapsearch m8 file; returns query, subject id, and e-value as a list"""
    m8_elts = m8_line.split('\t')
    return [m8_elts[0], m8_elts[1].replace(',',''), m8_elts[10]]
    # print m8_elts

def compare_evals(nr, mmetsp):
    """ compare evals between two searches """

    if float(mmetsp) < float(nr):
        return "True"
    else:
        return "False"

def compare_dicts(nr_dict, mmetsp_dict):
    """ Merge the two dictionaries."""


    EVAL_CUTOFF = -20

    match_counter = 0
    for query in nr_dict:
        if query in mmetsp_dict:
            nr_dict[query]["mmetsp_id"] = mmetsp_dict[query]["subject_id"]
            nr_dict[query]["mmetsp_e-value"] = mmetsp_dict[query]["e-value"]
            nr_dict[query]["mmetsp_match_best"] = compare_evals(nr_dict[query]["e-value"],nr_dict[query]["mmetsp_e-value"])
            del mmetsp_dict[query]
            match_counter += 1

        elif query not in mmetsp_dict:
            nr_dict[query]["mmetsp_id"] = "NA"
            nr_dict[query]["mmetsp_e-value"] = "NA"
            nr_dict[query]["mmetsp_match_best"] = "False"
    print match_counter, "double annotated contigs"
    return nr_dict

def write_out_csv(out_dict):
    """Write the results out in CSV format to the outfile"""

    for query in out_dict:
        out_elts = []
        out_elts.append(query)
        out_elts.append(out_dict[query]["subject_id"])
        out_elts.append(str(out_dict[query]["e-value"]))
        out_elts.append(out_dict[query]["mmetsp_id"])
        out_elts.append(str(out_dict[query]["mmetsp_e-value"]))
        out_elts.append(out_dict[query]["mmetsp_match_best"])

        outline = ",".join(out_elts)
        outline += "\n"
        outfile.write(outline)


print "Opening NR results file..."
nr_m8 = open(nr_m8_path, 'r')
print "Opening MMETSP results file..."
mmetsp_m8 = open(mmetsp_m8_path, 'r')
outfile = open(outfile_path, 'w')

print "Parsing NR m8 file..."
nr_dict = parse_m8(nr_m8)
nr_m8.close()
print "Parsing MMETSP m8 file..."
mmetsp_dict = parse_m8(mmetsp_m8)
mmetsp_m8.close()

print "NR matches:", len(nr_dict.keys())
print "MMETSP matches:", len(mmetsp_dict.keys())

print "Combining results..."
out_dict = compare_dicts(nr_dict, mmetsp_dict)

write_out_csv(out_dict)

outfile.close()
