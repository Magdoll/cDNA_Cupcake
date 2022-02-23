#!/usr/bin/env python
"""
A temporary CSV file for isoseq3 (v3.4+) dedup output

INPUT: dedup.fasta
OUTPUT: dedup.info.csv

2022 01 18 changes
    - allow for input/output file names while maintaining backward compat
    - allow for switched bc and umi positions in read header
"""
import os, re
import argparse
from Bio import SeqIO
from Bio.Seq import Seq

parser = argparse.ArgumentParser()
parser.add_argument('--input', '-i', default='dedup.fasta', help='Input dedup fasta file (default dedup.fasta)')
parser.add_argument('--output', '-o', default='dedup.info.csv', help='Output dedup csv (default dedup.info.csv)')
args = parser.parse_args()

rex = re.compile('(\S+) full_length_coverage=(\d+);length=(\d+);XM=(\S+);XC=(\S+)')
rex_switched = re.compile('(\S+) full_length_coverage=(\d+);length=(\d+);XC=(\S+);XM=(\S+)')
rex_umi_only = re.compile('(\S+) full_length_coverage=(\d+);length=(\d+);XM=(\S+)')

reader = SeqIO.parse(open(args.input),'fasta')
f = open(args.output, 'w')
f.write("id\tUMI\tUMIrev\tBC\tBCrev\tlength\tcount\n")
for r in reader:
    m = rex.match(r.description)
    m2 = rex_switched.match(r.description)
    if m is not None:
        _id, _count, _len, _umi, _bc = m.groups()
        f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(_id,_umi,Seq(_umi).reverse_complement(),_bc,Seq(_bc).reverse_complement(),_len,_count))
    elif not m and m2 is not None:
        _id, _count, _len, _bc, _umi = m2.groups()
        f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(_id,_umi,Seq(_umi).reverse_complement(),_bc,Seq(_bc).reverse_complement(),_len,_count))
    else: 
        m = rex_umi_only.match(r.description)
        _id, _count, _len, _umi = m.groups()
        _bc = 'NA'
        f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(_id,_umi,Seq(_umi).reverse_complement(),_bc,Seq(_bc).reverse_complement(),_len,_count))
f.close()
