#!/usr/bin/env python
__author__ = 'etseng@pacb.com'

"""
Helper script for calculating expected accuracy from FASTQ sequences.

ex:
i0_LQ_sample32f89e|c32/f1p0/929 isoform=c32;full_length_coverage=1;non_full_length_coverage=0;isoform_length=929;expected_accuracy=1.0
"""

import os, sys
from Bio import SeqIO

def phred_to_qv(phred):
    """Phred value to quality value."""
    return 10 ** -(phred / 10.0)

def calc_exp_acc(r, qv_trim_5, qv_trim_3):
    """
    :param r: SeqIO Fastq record
    :param qv_trim_5: 5' trimming
    :param qv_trim_3: 3' trimming
    """
    assert qv_trim_5 >= 0 and qv_trim_3 >= 0
    qv = r.letter_annotations['phred_quality']
    qv_len = len(qv)
    q = [phred_to_qv(x) for x in qv]
    if qv_trim_3 == 0:
        err_sum = sum(q[qv_trim_5:])
    else:
        err_sum = sum(q[qv_trim_5:-qv_trim_3])
    return 1.0 - (err_sum / float(qv_len))

def main(fastq_filename, output_filename, qv_trim_5, qv_trim_3):
    f = open(output_filename, 'w')
    for r in SeqIO.parse(open(fastq_filename), 'fastq'):
        exp_acc = calc_exp_acc(r, qv_trim_5, qv_trim_3)
        r.description += ";expected_accuracy={0:.3f}".format(exp_acc)
        SeqIO.write(r, f, 'fastq')
    f.close()



if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("fastq_filename",  help="FASTQ filename (ex: lq_isoforms.fastq")
    parser.add_argument("output_filename", help="Output FASTQ filename")
    parser.add_argument("--qv_trim_5", default=100, type=int, help="Ignore length on 5' for QV calculation (default: 100 bp)")
    parser.add_argument("--qv_trim_3", default=30, type=int, help="Ignore length on 3' for QV calculation (default: 30 bp)")

    args = parser.parse_args()
    main(args.fastq_filename, args.output_filename, args.qv_trim_5, args.qv_trim_3)
