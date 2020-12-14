#!/usr/bin/env python

import os, sys
from collections import defaultdict
from Bio import SeqIO
from cupcake.io import BioReaders

def type_fa_or_fq(file):
    file = file.upper()
    if file.endswith('.FA') or file.endswith('.FASTA'): return 'fasta'
    else: return 'fastq'

def summarize_GMAP_sam(input_fa_or_fq, input_sam):
    d = dict((r.id, len(r.seq)) for r in SeqIO.parse(open(input_fa_or_fq), type_fa_or_fq(input_fa_or_fq)))

    map_count = defaultdict(lambda: 0)
    for r in BioReaders.GMAPSAMReader(input_sam, True):
        map_count[r.qID] += 1
    multi = [x for x in map_count if map_count[x]>1]

    f = open(input_sam + '.summary.txt', 'w')
    f.write("id\tqLength\tqCoverage\tidentity\tnum_nonmatch\tnum_ins\tnum_del\tunique\n")
    for r in BioReaders.GMAPSAMReader(input_sam, True, query_len_dict=d):
        if r.sID == '*': continue
        if r.qID in multi: uni = 'N'
        else: uni = 'Y'
        f.write("{0}\t{1}\t{2:.4f}\t{3:.4f}\t\t{4}\t{5}\t{6}\t{7}\n".format(
            r.qID, d[r.qID], r.qCoverage, r.identity,
            r.num_nonmatches, r.num_ins, r.num_del,
            uni))
    f.close()

    print("Output written to: {0}".format(f.name), file=sys.stderr)

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser("Summarize GMAP SAM file in tab-delimited file.")
    parser.add_argument("input_fa_or_fq", help="Input fasta/fastq filename")
    parser.add_argument("sam_file", help="(GMAP) SAM filename")

    args = parser.parse_args()
    summarize_GMAP_sam(args.input_fa_or_fq, args.sam_file)
