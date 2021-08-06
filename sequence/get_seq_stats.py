#!/usr/bin/env python

__version__ = '1.0'

import os, sys
import gzip
import numpy as np
from Bio import SeqIO

def type_fa_or_fq(file):
    file = file.upper()
    if file.endswith('.FA') or file.endswith('.FASTA') or file.endswith('CLIPS'): return 'fasta'
    else: return 'fastq'

def get_seq_stats(file, binwidth):
    if file.endswith('.gz'):
        is_gzip = True
        filetype = type_fa_or_fq(file[:-3])
        reader = SeqIO.parse(gzip.open(file, 'rt'), filetype)
    else:
        is_gzip = False
        filetype = type_fa_or_fq(file)
        reader = SeqIO.parse(open(file), filetype)

    print("file type is:", filetype)

    f = open(file + '.seqlengths.txt', 'w')
    lens = []
    for r in reader:
        f.write(r.id + '\t' + str(len(r.seq)) + '\n')
        lens.append(len(r.seq))
    f.close()

    print("{0} sequences".format(len(lens)))
    print("min:", min(lens))
    print("max:", max(lens))
    print("avg:", sum(lens)*1./len(lens))

    # print by 1 kb bins
    print("Length Breakdown by kb range:")

    _max = (max(lens)//binwidth)+1
    bin = [0]*_max
    for x in lens: bin[x//binwidth] += 1

    for i in range(0, _max):
        if binwidth == 1000:
            print("{0}-{1} kb: {2}".format(i, i+1, bin[i]))
        else:
            print("{0}-{1}: {2}".format(i*binwidth, (i+1)*binwidth, bin[i]))

    print("5-95% percentile:", np.percentile(lens, 5), np.percentile(lens, 95))

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser("Summarize sequence lengths in fasta/fastq")
    parser.add_argument("filename", help="Input fasta/fastq filename, can be .gz")
    parser.add_argument("-b", "--binwidth", default=1000, type=int, help="Bin width, in bp (default: 1000 bp)")
    args = parser.parse_args()

    file = args.filename

    get_seq_stats(file, args.binwidth)

