#!/usr/bin/env python

__version__ = '1.0'

import os, sys
from Bio import SeqIO

def sort_by_len(input, reverse):
    if input.endswith('.fasta'):
        output = input[:-6] + '.sorted.fasta'
    elif input.endswith('.fa'):
        output = input[:-3] + '.sorted.fa'
    else:
        print >> sys.stderr, "Input must end with .fasta or .fa! Abort!"
        sys.exit(-1)
    
    with open(input) as f:
        d = SeqIO.to_dict(SeqIO.parse(f, 'fasta'))
    
    keys = d.keys()
    keys.sort(key=lambda x: len(d[x].seq),reverse=reverse)
    
    with open(output, 'w') as f:
        for k in keys: 
            f.write(">{0}\n{1}\n".format(k, d[k].seq))
    
    print >> sys.stderr, "Sorted output printed to", f.name

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser("Sort input fasta file")
    parser.add_argument("fasta_filename", help="Input fasta filename (must end with .fasta or .fa)")
    parser.add_argument("-r", "--reverse", default=False, action="store_true", help="Sort by decreasing length (default: off)")

    args = parser.parse_args()

    sort_by_len(args.fasta_filename, args.reverse)

