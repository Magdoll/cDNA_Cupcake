#!/usr/bin/env python
__author__ = 'etseng@pacb.com'

"""
Sequence headers must have:

@i0_LQ_sample32f89e|c31/f12p0/433 isoform=c31;full_length_coverage=12;non_full_length_coverage=0;isoform_length=433;expected_accuracy=0.300

"""
import os, sys
from Bio import SeqIO


def func(r):
    fl_count = None
    exp_acc = None
    for x in r.description.split(';'):
        if x.find('=')==-1: continue
        a,b=x.split('=')
        if a=='full_length_coverage':
            fl_count = int(b)
        elif a=='expected_accuracy':
            exp_acc = float(b)
    return fl_count, exp_acc


def main(fastq_filename, output_filename, min_fl_count, min_exp_acc, is_flnc):
    f = open(output_filename, 'w')
    for r in SeqIO.parse(open(fastq_filename), 'fastq'):
        fl_count, exp_acc = func(r)
        if not is_flnc and fl_count is None:
            print("Sequence header does not include field `full_length_coverage=`. Abort!", file=sys.stderr)
            sys.exit(-1)
        if exp_acc is None:
            print("Sequence header does not include field `expected_accuracy=`. Please run calc_expected_accuracy_from_fastq.py script first!", file=sys.stderr)
            sys.exit(-1)
        if (is_flnc or fl_count >= min_fl_count) and exp_acc >= min_exp_acc:
            print("Including {0} into output file (FL: {1}, acc: {2}).".format(r.id, fl_count, exp_acc), file=sys.stderr)
            SeqIO.write(r, f, 'fastq')
    f.close()

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("fastq_filename",  help="LQ FASTQ filename (ex: lq_isoforms.fastq")
    parser.add_argument("output_filename", help="Output FASTQ filename")
    parser.add_argument("--min_fl_count", default=2, type=int, help="Minimum FL count (default: 2).")
    parser.add_argument("--min_exp_acc", type=float, default=0.99, help="Minimum predicted accuracy (default: 0.99).")
    parser.add_argument("--is_flnc", action="store_true", default=False, help="Input FASTQ is FLNC, not LQ")

    args = parser.parse_args()
    assert 0 <= args.min_exp_acc <= 1
    main(args.fastq_filename, args.output_filename, args.min_fl_count, args.min_exp_acc, args.is_flnc)
