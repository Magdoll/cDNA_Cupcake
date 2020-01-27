#!/usr/bin/env python
import os, sys
from cupcake.io import GFF

def get_gff_from_list(gff_filename, listfile, partial_ok=False):
    seqs = [line.strip() for line in open(listfile)]
    for r in GFF.collapseGFFReader(gff_filename):
        if r.seqid in seqs or r.seqid.split('|')[0] in seqs or (partial_ok and any(r.seqid.startswith(x) for x in seqs)):
            GFF.write_collapseGFF_format(sys.stdout, r)

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser("Get records from a GFF file from a list")
    parser.add_argument("gff_filename", help="Input gff filename to extract sequences from")
    parser.add_argument("list_filename", help="List of sequence IDs to extract")
    parser.add_argument("--partial", action="store_true", default=False, help="OK if seq IDs only match the beginning")

    args = parser.parse_args()

    get_gff_from_list(args.gff_filename, args.list_filename, args.partial)
