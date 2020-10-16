#!/usr/bin/env python
"""
A temporary CSV file for isoseq3 (v3.4+) dedup output

INPUT: dedup.fasta
OUTPUT: dedup.info.csv
"""
import os, re
from Bio import SeqIO
from Bio.Seq import Seq

rex = re.compile('(\S+) full_length_coverage=(\d+);length=(\d+);XM=(\S+);XC=(\S+)')
rex_umi_only = re.compile('(\S+) full_length_coverage=(\d+);length=(\d+);XM=(\S+)')

reader = SeqIO.parse(open('dedup.fasta'),'fasta')
f = open('dedup.info.csv', 'w')
f.write("id\tUMI\tUMIrev\tBC\tBCrev\tlength\tcount\n")
for r in reader:
    m = rex.match(r.description)
    if m is not None:
        _id, _count, _len, _umi, _bc = m.groups()
        f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(_id,_umi,Seq(_umi).reverse_complement(),_bc,Seq(_bc).reverse_complement(),_len,_count))
    else: 
        m = rex_umi_only.match(r.description)
        _id, _count, _len, _umi = m.groups()
        _bc = 'NA'
        f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(_id,_umi,Seq(_umi).reverse_complement(),_bc,Seq(_bc).reverse_complement(),_len,_count))
f.close()
