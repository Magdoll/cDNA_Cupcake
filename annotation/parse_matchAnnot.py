#!/usr/bin/env python
import os, sys
from collections import defaultdict
from Bio import SeqIO

def type_fa_or_fq(file):
    file = file.upper()
    if file.endswith('.FA') or file.endswith('.FASTA'): return 'fasta'
    else: return 'fastq'

def parse_matchAnnot(fa_or_fq, filename, not_pbid=False, parse_FL_coverage=False, fusion_mode=False):

    pbids = []
    fl_cov = {} # only used if parse_FL_coverage is True
    for r in SeqIO.parse(open(fa_or_fq), type_fa_or_fq(fa_or_fq)):
        _id = r.id if not_pbid else r.id.split('|')[0]
        pbids.append(_id)
        if parse_FL_coverage:
            try:
                cov = int(r.description.split('full_length_coverage=')[1].split(';')[0])
                fl_cov[_id] = cov
            except:
                print >> sys.stderr, "WARNING: Unable to extract `full_length_coverage=` from {0}. Mark as NA.".format(r.description)
                fl_cov[_id] = 'NA'

    if not fusion_mode:
        match = defaultdict(lambda: (None,None,0)) # ex: PB.1.1 -> (NOC2L, NOC2L-001, 5)
    else:
        match = defaultdict(lambda: [])

    for line in open(filename):
        if line.startswith('result:'):
            raw = line.strip().split()
            if len(raw) < 7: continue
            pbid = raw[1] if not_pbid else raw[1].split('|')[0]
            gene = raw[2]
            isoform = raw[3]
            score = int(raw[7])
            if not fusion_mode: # keep the highest score
                if score > match[pbid][1]: match[pbid] = (gene, isoform, score)
            else:
                match[pbid].append((gene, isoform, score))

    f = open(filename + '.parsed.txt', 'w')
    f.write("pbid\tpbgene\trefisoform\trefgene\tscore")
    if parse_FL_coverage: f.write("\tcount_fl")
    f.write("\n")
    for pbid in pbids:
        if not_pbid:
            pbpre = pbid
        else:
            pbpre = pbid.split('.')[1]
        _cov_text = "\t{0}".format(fl_cov[pbid]) if parse_FL_coverage else ''
        if pbid not in match:
            f.write("{0}\t{1}\tNA\tNA\tNA{2}\n".format(pbid, pbpre, _cov_text))
        else:
            if not fusion_mode:
                gene, isoform, score = match[pbid]
                if gene is None: f.write("{0}\t{1}\tNA\tNA\tNA{2}\n".format(pbid, pbpre, _cov_text))
                else:
                    f.write("{0}\t{1}\t{2}\t{3}\t{4}{5}\n".format(pbid, pbpre, isoform, gene, score, _cov_text))
            else:
                results = match[pbid]
                if len(results) == 0: f.write("{0}\t{1}\tNA\tNA\tNA{2}\n".format(pbid, pbpre, _cov_text))
                else:
                    genes = ",".join(x[0] for x in results)
                    isoforms = ",".join(x[1] for x in results)
                    scores = ",".join(str(x[2]) for x in results)  
                    f.write("{0}\t{1}\t{2}\t{3}\t{4}{5}\n".format(pbid, pbpre, isoforms, genes, scores, _cov_text))
    f.close()
    print >> sys.stderr, "Output written to:", f.name

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser("Parse MatchAnnot result")
    parser.add_argument("fa_or_fq", help="Fasta/Fastq filename used to create the SAM file for matchAnnot")
    parser.add_argument("match_filename", help="MatchAnnot filename")
    parser.add_argument("--not_pbid", action="store_true", default=False, help="Turn this on if not sequence ID is not PB.X.Y (default: off)")
    parser.add_argument("--parse_FL_coverage", action="store_true", default=False, help="Parse `full_length_coverage=` from sequence ID.")
    parser.add_argument("--fusion_mode", action="store_true", default=False, help="Input is fusion SAM mapping (default: off)")

    args = parser.parse_args()
    parse_matchAnnot(args.fa_or_fq, args.match_filename, args.not_pbid, args.parse_FL_coverage, args.fusion_mode)

