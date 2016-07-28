#!/usr/bin/env python

__version__ = '1.0'

import os, sys
import numpy as np
from Bio import SeqIO

def type_fa_or_fq(file):
    file = file.upper()
    if file.endswith('.FA') or file.endswith('.FASTA'): return 'fasta'
    else: return 'fastq'

def get_seq_stats(file):
    print "file type is:", type_fa_or_fq(file)

    f = open(file + '.seqlengths.txt', 'w')
    lens = []
    for r in SeqIO.parse(open(file), type_fa_or_fq(file)):
        f.write(r.id + '\t' + str(len(r.seq)) + '\n')
        lens.append(len(r.seq))
    f.close()


    print "{0} sequences".format(len(lens))
    print "min:", min(lens)
    print "max:", max(lens)
    print "avg:", sum(lens)*1./len(lens)

    # print by 1 kb bins
    print "Length Breakdown by kb range:"

    _max = max(lens)/1000+1
    bin = [0]*_max
    for x in lens: bin[x/1000] += 1

    for i in xrange(0, _max):
        print "{0}-{1} kb: {2}".format(i, i+1, bin[i])

    print "5-95% percentile:", np.percentile(lens, 5), np.percentile(lens, 95)

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser("Summarize sequence lengths in fasta/fastq")
    parser.add_argument("filename", help="Input fasta/fastq filename")
    args = parser.parse_args()

    file = args.filename

    get_seq_stats(file)

