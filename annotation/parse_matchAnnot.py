#!/usr/bin/env python
import os, sys
from collections import defaultdict
from Bio import SeqIO

def parse_matchAnnot(fq_filename, filename, not_pbid=False):
    pbids = [r.id.split('|')[0] for r in SeqIO.parse(open(fq_filename), 'fastq')]
    match = defaultdict(lambda: (None,None,0)) # ex: PB.1.1 -> (NOC2L, NOC2L-001, 5)

    for line in open(filename):
        if line.startswith('result:'):
            raw = line.strip().split()
            if len(raw) < 7: continue
            pbid = raw[1].split('|')[0]
            gene = raw[2]
            isoform = raw[3]
            score = int(raw[7])
            if score > match[pbid][1]: match[pbid] = (gene, isoform, score)
    
    f = open(filename + '.parsed.txt', 'w')
    f.write("pbid\tpbgene\trefisoform\trefgene\tscore\n")
    for pbid in pbids:
        if not_pbid:
            pbpre = pbid
        else:
            pbpre = pbid.split('.')[1]
        if pbid not in match: f.write("{0}\t{1}\tNA\tNA\tNA\n".format(pbid, pbpre))
        else:
            gene, isoform, score = match[pbid]
            if gene is None: f.write("{0}\t{1}\tNA\tNA\tNA\n".format(pbid, pbpre))
            else:
                f.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(pbid, pbpre, isoform, gene, score))
    f.close()
    print >> sys.stderr, "Output written to:", f.name

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser("Parse MatchAnnot result")
    parser.add_argument("fastq_filename", help="Fastq filename used to create the SAM file for matchAnnot")
    parser.add_argument("match_filename", help="MatchAnnot filename")
    parser.add_argument("--not_pbid", action="store_true", default=False, help="Turn this on if not sequence ID is not PB.X.Y (default: off)")

    args = parser.parse_args()
    parse_matchAnnot(args.fastq_filename, args.match_filename, args.not_pbid)

