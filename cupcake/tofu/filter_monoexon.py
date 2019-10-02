#!/usr/bin/env python

__author__ = 'etseng@pacb.com'

"""
Filter away transcripts that are mono-exonic.

Required input: <input_prefix>.gff, .rep.fq, .abundance.txt

Example:
    filter_monoexon.py test.collapsed test.collapsed.nomono

"""

import os, sys
from collections import defaultdict
from csv import DictReader, DictWriter

from Bio import SeqIO

from cupcake.io import GFF


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
        print("File {0} does not exist. Abort!".format(count_filename), file=sys.stderr)
        sys.exit(-1)
    if not os.path.exists(gff_filename):
        print("File {0} does not exist. Abort!".format(gff_filename), file=sys.stderr)
        sys.exit(-1)
    if not os.path.exists(rep_filename):
        print("File {0} does not exist. Abort!".format(rep_filename), file=sys.stderr)
        sys.exit(-1)

    pbids1 = set([r.id for r in SeqIO.parse(open(rep_filename),'fastq')])
    pbids2 = set([r.seqid for r in GFF.collapseGFFReader(gff_filename)])
    pbids3 = set(read_count_file(count_filename)[0].keys())

    if len(pbids1)!=len(pbids2) or len(pbids2)!=len(pbids3) or len(pbids1)!=len(pbids3):
        print("The number of PBID records in the files disagree! Sanity check failed.", file=sys.stderr)
        print("# of PBIDs in {0}: {1}".format(rep_filename, len(pbids1)), file=sys.stderr)
        print("# of PBIDs in {0}: {1}".format(gff_filename, len(pbids2)), file=sys.stderr)
        print("# of PBIDs in {0}: {1}".format(count_filename, len(pbids3)), file=sys.stderr)
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


def main():
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("input_prefix", help="Input prefix (ex: test.collapsed.min_fl_2)")

    args = parser.parse_args()
    output_prefix = args.input_prefix + '.nomono'

    count_filename, gff_filename, rep_filename = sanity_check_collapse_input(args.input_prefix)

    good = []
    f = open(output_prefix + '.gff', 'w')
    reader = GFF.collapseGFFReader(gff_filename)
    for r in reader:
        assert r.seqid.startswith('PB.')
        if len(r.ref_exons) > 1:
            good.append(r.seqid)
            GFF.write_collapseGFF_format(f, r)

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

    print("Output written to:", output_prefix + '.gff', file=sys.stderr)
    print("Output written to:", output_prefix + '.rep.fq', file=sys.stderr)
    print("Output written to:", output_prefix + '.gff', file=sys.stderr)

if __name__ == "__main__":
    main()
