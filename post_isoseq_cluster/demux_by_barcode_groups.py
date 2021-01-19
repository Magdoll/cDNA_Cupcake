#!/usr/bin/env python
__author__="etseng@pacb.com"

"""
Given a pooled input GFF + demux CSV file, write out per-{barcode group} GFFs
If input fasta/fastq is given, optionally also output per-{barcode group} FASTA/FASTQ
"""
import os, re, sys
import cupcake.io.GFF as GFF
from csv import DictReader
from collections import defaultdict
from Bio import SeqIO

rex_pbid = re.compile('(PB.\d+.\d+)(|\S+)')

def get_type_fafq(in_filename):
    in_filename = in_filename.upper()
    if in_filename.endswith('.FA') or in_filename.endswith('FASTA'): return 'fasta'
    elif in_filename.endswith('.FQ') or in_filename.endswith('FASTQ'): return 'fastq'
    else:
        raise Exception("Unrecognized file suffix .{0}! Must end with .fasta or .fastq!".format(in_filename[in_filename.find('.'):]))

def regroup_gff(pooled_gff, demux_count_file, output_prefix, out_group_dict, in_fafq=None):
    """
    :param pooled_sam: SAM file
    :param demux_count_file: comma-delimited per-barcode count file
    :param output_prefix: output prefix for GFF
    :param out_group_dict: dict of barcode name --> group to be long in  (ex: {'EM1':'EM', 'EM2':'EM'})
    :param in_fafq: optional fasta/fastq that was input to SAM
    """
    if in_fafq is not None: type_fafq = get_type_fafq(in_fafq)
    in_tissue = defaultdict(lambda: set()) # pbid --> list of tissue it is in (EM, END, R)

    for r in DictReader(open(demux_count_file),delimiter=','):
        for k,v in r.items():
            if k!='id' and int(v) > 0: in_tissue[r['id']].add(k)

#    in_tissue = dict(in_tissue)

    handles = {}
    handles_fafq = {}
    for g in out_group_dict.values():
        handles[g] = open("{o}_{g}_only.gff".format(o=output_prefix, g=g), 'w')
        if in_fafq is not None: handles_fafq[g] = open("{o}_{g}_only.{t}".format(o=output_prefix, g=g, t=type_fafq), 'w')

    if in_fafq is not None:
        fafq_dict = SeqIO.to_dict(SeqIO.parse(open(in_fafq), type_fafq))
        fafq_dict_keys = list(fafq_dict.keys())
        for k in fafq_dict_keys:
            m = rex_pbid.match(k)
            if m is not None: fafq_dict[m.group(1)] = fafq_dict[k]
    reader = GFF.collapseGFFReader(pooled_gff)
    for r in reader:
        groups_to_write_in = set()
        pbid = r.seqid
        if pbid not in in_tissue:
            print("WARNING: {0} does not belong to any group indicated by outgroup_dict".format(pbid), file=sys.stderr)
        for tissue in in_tissue[pbid]:
            groups_to_write_in.add(out_group_dict[tissue])

        for g in groups_to_write_in:
            GFF.write_collapseGFF_format(handles[g], r)
            if in_fafq is not None:
                SeqIO.write(fafq_dict[pbid], handles_fafq[g], type_fafq)

if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("--pooled_fastx", default=None, help="Pooled FASTA/FASTQ (optional, if given, will also output demux fa/fq)")
    parser.add_argument("pooled_gff", help="Pooled GFF file")
    parser.add_argument("demux_count_file", help="Demux count file")
    parser.add_argument("output_prefix", help="Output prefix for GFF outputs")
    parser.add_argument("outgroup_dict", help="Tuples indicating barcode grouping")

    args = parser.parse_args()
    tmp = eval(args.outgroup_dict)
    out_group_dict = dict([tmp]) if len(tmp)==1 else dict(tmp)
    regroup_gff(args.pooled_gff, args.demux_count_file, args.output_prefix, out_group_dict, args.pooled_fastx)
