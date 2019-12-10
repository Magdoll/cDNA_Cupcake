#!/usr/bin/env python
__author__="etseng@pacb.com"

import os, re, sys
from Bio import SeqIO
from csv import DictReader, DictWriter
from collections import Counter

FIELDNAMES = ['pbid', 'pbgene', 'length', 'refisoform', 'refgene', 'fl_count']
pbid_rex = re.compile('PB.(\d+).\d+')

def demux_for_subsamping(class_filename, fasta_filename, demux_count_file, output_prefix, out_group_dict, ignore_novel):
    # read SQANTI classification to get known gene/transcript name
    d = {} # pbid --> record
    for r in DictReader(open(class_filename),delimiter='\t'):
        d[r['isoform']] = r

    # get read lengths
    lens = {} # pbid -> length
    for r in SeqIO.parse(open(fasta_filename),'fasta'):
        lens[r.id] = len(r.seq)

    writers = {}
    handles = {}
    out_groups = set(out_group_dict.values())
    for g in out_groups:
        handles[g] = open("{o}_{g}_only.{i}.for_subsampling.txt".format(\
            o=output_prefix, g=g, i="ignore_novel" if ignore_novel else "use_novel"), 'w')
        writers[g] = DictWriter(handles[g], FIELDNAMES, delimiter='\t')
        writers[g].writeheader()


    reader = DictReader(open(demux_count_file),delimiter=',')
    for r in reader:
        if r['id'] not in d:
            print("WARNING: skipping {0} because not in {1}".format(\
                r['id'], class_filename), file=sys.stderr)
            continue

        m = pbid_rex.match(r['id'])
        if m is None:
            print("ERROR: unable to parse ID {0}. Expected format PB.X.Y!".format(r['id']), file=sys.stderr)
            sys.exit(-1)

        newrec = {'pbid': r['id'], 'pbgene': m.group(1), 'length': lens[r['id']]}

        gene = d[r['id']]['associated_gene']
        trans = d[r['id']]['associated_transcript']
        if gene.startswith('novel') and ignore_novel: gene = 'NA'
        if trans.startswith('novel'):
            if ignore_novel: trans = 'NA'
            else: trans += r['id'] # add an unique identified to make this "novel" refgene unique
        newrec['refgene'] = gene
        newrec['refisoform'] = trans

        group_counts = Counter()
        for b, g in out_group_dict.items():
            group_counts[g] += int(r[b])

        for g in out_groups:
            newrec['fl_count'] = group_counts[g]
            writers[g].writerow(newrec)

    for h in handles.values():
        h.close()

if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("class_filename", help="SQANTI classification file")
    parser.add_argument("fasta_filename", help="FASTA filename")
    parser.add_argument("demux_count_file", help="Demux count file")
    parser.add_argument("output_prefix", help="Output prefix for GFF outputs")
    parser.add_argument("outgroup_dict", help="Tuples indicating barcode grouping")
    parser.add_argument("--ignore_novel", action="store_true", default=False, help="Ignore novel genes/transcripts (default: off)")


    args = parser.parse_args()
    out_group_dict = dict(eval(args.outgroup_dict))
    demux_for_subsamping(args.class_filename, args.fasta_filename, args.demux_count_file, args.output_prefix, out_group_dict, args.ignore_novel)
