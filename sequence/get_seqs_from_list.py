#!/usr/bin/env python
import os, sys
from Bio import SeqIO

def get_seqs_from_list(fastafile, listfile):
    seqs = [line.strip() for line in open(listfile)]
    for r in SeqIO.parse(open(fastafile), 'fasta'):
        if r.id in seqs or r.id.split('|')[0] in seqs or any(r.id.startswith(x) for x in seqs):
            print ">" + r.id
            print r.seq

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser("Get sequences from a fasta file from a list")
    parser.add_argument("fasta_filename", help="Input fasta filename to extract sequences from")
    parser.add_argument("list_filename", help="List of sequence IDs to extract")

    args = parser.parse_args()

    get_seqs_from_list(args.fasta_filename, args.list_filename)
