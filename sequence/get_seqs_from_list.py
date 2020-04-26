#!/usr/bin/env python
import os, sys
from Bio import SeqIO

def get_seqs_from_list(fafq, listfile, partial_ok=False, exclude=False):
    f = sys.stdout
    type = 'fastq' if fafq.upper().endswith('.FQ') or fafq.upper().endswith('.FASTQ') else 'fasta'
    seqs = set([line.strip() for line in open(listfile)])
    for r in SeqIO.parse(open(fafq), type):
        id_seen = r.id in seqs or r.id.split('|')[0] in seqs or (partial_ok and any(r.id.startswith(x) for x in seqs))
        if id_seen ^ exclude:
            SeqIO.write(r, f, type)

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser("Get sequences from a fasta/fastq file from a list")
    parser.add_argument("fasta_filename", help="Input fasta/fastq filename to extract sequences from")
    parser.add_argument("list_filename", help="List of sequence IDs to extract")
    parser.add_argument("--partial", action="store_true", default=False, help="OK if seq IDs only match the beginning")
    parser.add_argument("--exclude", action="store_true", default=False, help="Output sequences NOT in the list, default OFF")

    args = parser.parse_args()

    get_seqs_from_list(args.fasta_filename, args.list_filename, args.partial, args.exclude)
