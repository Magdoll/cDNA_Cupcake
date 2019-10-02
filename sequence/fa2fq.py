#!/usr/bin/env python

__version__ = '1.0'

from Bio import SeqIO
import os, sys

def main():
    from argparse import ArgumentParser
    parser = ArgumentParser("Convert fasta to fastq")
    parser.add_argument("fasta_filename", help="input fasta (must end with .fasta or .fa)")
    args = parser.parse_args()

    input = args.fasta_filename

    fa2fq(input)

def fa2fq(input):
    try:
        assert input.lower().endswith('.fasta') or input.lower().endswith('.fa')
    except AssertionError:
        print("Input {0} does not end with .fasta or .fa! Abort".format(input), file=sys.stderr)
        sys.exit(-1)
    output = input[:input.rfind('.')] + '.fastq'

    f = open(output, 'w')
    for r in SeqIO.parse(open(input),'fasta'):
        r.letter_annotations['phred_quality'] = [60]*len(r.seq)
        SeqIO.write(r, f, 'fastq')
    f.close()

    print("Output written to", f.name, file=sys.stderr)
    return f.name

if __name__ == "__main__":
    main()
