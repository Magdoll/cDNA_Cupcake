#!/usr/bin/env python

__version__ = '1.0'


import os, sys, random
from collections import defaultdict
from Bio import SeqIO

def type_fa_or_fq(file):
    file = file.upper()
    if file.endswith('.FA') or file.endswith('.FASTA'): return 'fasta'
    else: return 'fastq'

def sep_by_primer(file, output_prefix, sample_size):
    filetype = type_fa_or_fq(file)

    ids = [r.id for r in SeqIO.parse(open(file), filetype)]

    n = len(ids)
    if sample_size > n:
        print("WARNING: {0} contains only {1} sequences but subsample size at {2}! Simply output whole file.".format(\
            file, n, sample_size), file=sys.stderr)

    chosen_ids = random.sample(ids, min(n, sample_size))

    with open(output_prefix+'.'+'random'+str(sample_size)+'.'+filetype, 'w') as f:
        for r in SeqIO.parse(open(file), filetype):
            if r.id in chosen_ids:
                SeqIO.write(r, f, filetype)

        print("Randomly selected sequences written to {0}.".format(f.name), file=sys.stderr)


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser("Randomly select N sequences from fasta/fastq files")
    parser.add_argument("filename", help="Input fasta/fastq filename")
    parser.add_argument("output_prefix", help="Output file prefix")
    parser.add_argument("sample_size", type=int, help="Subsample size")
    args = parser.parse_args()

    sep_by_primer(args.filename, args.output_prefix, args.sample_size)


