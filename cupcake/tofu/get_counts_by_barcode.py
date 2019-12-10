#!/usr/bin/env python
__author__ = 'etseng@pacb.com'

"""
Used for getting per-barcode Fl count information after IsoSeq cluster is run and HQ is mapped and collapsed.

Required Input:

isoseq_classify.primer_info.csv -- CSV output from classify
*.collapsed.group.txt, *.collapsed.rep.fq --- to get the PB IDs and group info

Optional input: primer names in a text file where each line is:

<primer_index>, <primer_name>

"""

import os, sys
from csv import DictReader
from collections import defaultdict
from collections import Counter


def read_classify_csv(csv_filename):
    """
    :param csv_filename: CSV file from IsoSeq classify, must have "id" and "primer" field
    :return: primer_ranges (list of indices), primer_info which is dict of read_id --> primer (all strings)
    """
    primer_info = {}
    for r in DictReader(open(csv_filename),delimiter=','):
        if r['primer']!='NA': primer_info[r['id']] = r['primer']
    primer_ranges = list(set(map(int, list(primer_info.values()))))
    primer_ranges.sort()
    primer_ranges = list(map(str, primer_ranges))
    return primer_ranges, primer_info


#primer_names = ['B6-GV', 'B6-MII', 'B6-1C', 'B6-2C', 'B6-8C', 'B6-BI', 'DBA-GV', 'DBA-MII', 'DBA-1C', 'DBA-2C', 'DBA-8C', 'DBA-BI']

def get_fl_count_by_barcode(collapse_prefix, classify_csv, cluster_csv, primer_names=None):

    primer_ranges, primer_info = read_classify_csv(classify_csv)
    if primer_names is None:
        primer_names = dict((int(x), x) for x in primer_ranges)

    cluster_info = defaultdict(lambda: [])
    for r in DictReader(open(cluster_csv),delimiter=','):
        cluster_info[r['cluster_id']].append(r)

    group_filename = collapse_prefix + '.group.txt'
    print("Reading {0}....".format(group_filename), file=sys.stderr)

    f = open(collapse_prefix + '.fl_count_by_barcode.txt', 'w')
    f.write("pbid")
    for p in primer_ranges: f.write('\t' + str(primer_names[int(p)]))
    f.write('\n')
    for line in open(group_filename):
        pbid, members = line.strip().split('\t')
        tally = Counter()
        for m in members.split(','):
            cid = m.split('/')[0]
            for r in cluster_info[cid]:
                if r['read_type'] == 'FL':
                    tally[primer_info[r['read_id']]] += 1
        f.write(pbid)
        for p in primer_ranges: f.write('\t' + str(tally[p]))
        f.write('\n')
    f.close()
    print("Output written to: {0}.".format(f.name), file=sys.stderr)


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("collapse_prefix", help="Collapse prefix (ex: hq_isoforms.fastq.no5merge.collapsed)")
    parser.add_argument("classify_csv", help="Classify output CSV (ex: classify.primer_info.csv)")
    parser.add_argument("cluster_csv", help="Cluster output CSV (ex: cluster_report.csv)")

    args = parser.parse_args()

    get_fl_count_by_barcode(args.collapse_prefix, args.classify_csv, args.cluster_csv)
