#!/usr/bin/env python

__author__ = 'etseng@pacb.com'

"""
Filter away 5' degraded products (isoforms that are shorter isoforms of longer ones)
3' differences are always honored.

If input is:

Isoform A: exon 1, 2, 3, 4
Isoform B: exon 3, 4
Isoform C: exon 3, 4, 5

Then Isoform A and C are preserved. Isoform B (degraded form of A) is filtered out.

Required input: <input_prefix>.gff, .rep.fq, .abundance.txt

Example:
    filter_away_subset.py test.collapsed test.collapsed.filtered

"""

import os, sys
from collections import defaultdict
from csv import DictReader, DictWriter

from Bio import SeqIO

from cupcake.io import GFF
from cupcake.tofu import compare_junctions

def sanity_check_collapse_input(input_prefix):
    """
    Check that
    1. the count, gff, rep files exist
    2. the number of records agree among the three
    """
    group_filename = input_prefix + '.group.txt'
    count_filename = input_prefix + '.abundance.txt'
    gff_filename = input_prefix + '.gff'
    rep_filename = input_prefix + '.rep.fq'
    if not os.path.exists(count_filename):
        print >> sys.stderr, "File {0} does not exist. Abort!".format(count_filename)
        sys.exit(-1)
    if not os.path.exists(gff_filename):
        print >> sys.stderr, "File {0} does not exist. Abort!".format(gff_filename)
        sys.exit(-1)
    if not os.path.exists(rep_filename):
        print >> sys.stderr, "File {0} does not exist. Abort!".format(rep_filename)
        sys.exit(-1)

    pbids1 = set([r.id for r in SeqIO.parse(open(rep_filename),'fastq')])
    pbids2 = set([r.seqid for r in GFF.collapseGFFReader(gff_filename)])
    pbids3 = set(read_count_file(count_filename)[0].keys())

    if len(pbids1)!=len(pbids2) or len(pbids2)!=len(pbids3) or len(pbids1)!=len(pbids3):
        print >> sys.stderr, "The number of PBID records in the files disagree! Sanity check failed."
        print >> sys.stderr, "# of PBIDs in {0}: {1}".format(rep_filename, len(pbids1))
        print >> sys.stderr, "# of PBIDs in {0}: {1}".format(gff_filename, len(pbids2))
        print >> sys.stderr, "# of PBIDs in {0}: {1}".format(count_filename, len(pbids3))
        sys.exit(-1)

    return count_filename, gff_filename, rep_filename


def read_count_file(count_filename):
    f = open(count_filename)
    count_header = ''
    while True:
        cur_pos = f.tell()
        line = f.readline()
        if not line.startswith('#'):
            f.seek(cur_pos)
            break
        else:
            count_header += line
    d = dict((r['pbid'], r) for r in DictReader(f, delimiter='\t'))
    f.close()
    return d, count_header


def can_merge(m, r1, r2, internal_fuzzy_max_dist):
    if m == 'subset':
        r1, r2 = r2, r1 #  rotate so r1 is always the longer one
    if m == 'super' or m == 'subset':
        n2 = len(r2.ref_exons)
        if r1.strand == '+':
            # if r2 is monoexonic, it can start after the last exon of r1's last exon
            # if r2 is multiexonic, the last start must be pretty close (fuzzy allowed)
            if n2==1: # r2 is mono-exonic
                return r1.ref_exons[-1].start - r2.ref_exons[-1].start <= internal_fuzzy_max_dist 
            else: return abs(r1.ref_exons[-1].start - r2.ref_exons[-1].start) <= internal_fuzzy_max_dist and \
                    r1.ref_exons[-n2].start <= r2.ref_exons[0].start < r1.ref_exons[-n2].end
        else:
            if n2==1: return r1.ref_exons[0].end - r2.ref_exons[0].end >= -internal_fuzzy_max_dist
            else: return abs(r1.ref_exons[0].end - r2.ref_exons[0].end) <= internal_fuzzy_max_dist and \
                    r1.ref_exons[n2-1].start <= r2.ref_exons[-1].end < r1.ref_exons[n2].end

def filter_out_subsets(recs, internal_fuzzy_max_dist):
    # recs must be sorted by start becuz that's the order they are written
    i = 0
    while i < len(recs)-1:
        no_change = True
        j = i + 1
        while j < len(recs):
            if recs[j].start > recs[i].end: 
                break
            recs[i].segments = recs[i].ref_exons
            recs[j].segments = recs[j].ref_exons
            m = compare_junctions.compare_junctions(recs[i], recs[j], internal_fuzzy_max_dist)
            if can_merge(m, recs[i], recs[j], internal_fuzzy_max_dist):
                if m == 'super': # pop recs[j] 
                    recs.pop(j)
                else:
                    recs.pop(i)
                    no_change = False
            else:
                j += 1
        if no_change: i += 1


def main():
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("input_prefix", help="Input prefix (ex: test.collapsed.min_fl_2)")
    parser.add_argument("--fuzzy_junction", type=int, default=5, help="Fuzzy junction max dist (default: 5bp)")

    args = parser.parse_args()
    output_prefix = args.input_prefix + '.filtered'

    count_filename, gff_filename, rep_filename = sanity_check_collapse_input(args.input_prefix)

    recs = defaultdict(lambda: [])
    reader = GFF.collapseGFFReader(gff_filename)
    for r in reader:
        assert r.seqid.startswith('PB.')
        recs[int(r.seqid.split('.')[1])].append(r)

    good = []
    f = open(output_prefix + '.gff', 'w')
    keys = recs.keys()
    keys.sort()
    for k in recs:
        xxx = recs[k]
        filter_out_subsets(xxx, args.fuzzy_junction)
        for r in xxx:
            GFF.write_collapseGFF_format(f, r)
            good.append(r.seqid)
    f.close()

    # read abundance first
    d, count_header = read_count_file(count_filename)

    # write output rep.fq
    f = open(output_prefix + '.rep.fq', 'w')
    for r in SeqIO.parse(open(rep_filename), 'fastq'):
        if r.name.split('|')[0] in good:
            SeqIO.write(r, f, 'fastq')
    f.close()

    # write output to .abundance.txt
    f = open(output_prefix + '.abundance.txt', 'w')
    f.write(count_header)
    writer = DictWriter(f, fieldnames=['pbid','count_fl','count_nfl','count_nfl_amb','norm_fl','norm_nfl','norm_nfl_amb'], \
                        delimiter='\t', lineterminator='\n')
    writer.writeheader()
    for k in good:
        r = d[k]
        writer.writerow(r)
    f.close()

    print >> sys.stderr, "Output written to:", output_prefix + '.gff'
    print >> sys.stderr, "Output written to:", output_prefix + '.rep.fq'
    print >> sys.stderr, "Output written to:", output_prefix + '.gff'

if __name__ == "__main__":
    main()
