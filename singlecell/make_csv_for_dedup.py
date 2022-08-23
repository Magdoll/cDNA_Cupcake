#!/usr/bin/env python
"""
A temporary CSV file for isoseq3 (v3.4+) dedup output

INPUT: dedup.fasta
OUTPUT: dedup.info.csv
"""
import os, re
from Bio import SeqIO
from Bio.Seq import Seq
import pysam

rex = re.compile('(\S+) full_length_coverage=(\d+);length=(\d+);XM=(\S+);XC=(\S+)')
rex_umi_only = re.compile('(\S+) full_length_coverage=(\d+);length=(\d+);XM=(\S+)')

if os.path.exists('dedup.bam'):
    _type = 'BAM'
    reader = pysam.AlignmentFile(open('dedup.bam'),'rb',check_sq=False)
elif os.path.exists('dedup.fasta'):
    _type = 'FASTA'
    reader = SeqIO.parse(open('dedup.fasta'),'fasta')
else:
    print("Cannot find dedup.bam or dedup.fasta. Abort!")
    sys.exit(-1)
f = open('dedup.info.csv', 'w')
f.write("id\tUMI\tUMIrev\tBC\tBCrev\tlength\tcount\n")
for r in reader:
    if _type == 'FASTA':
         m = rex.match(r.description)
         _id, _count, _len, _umi, _bc = m.groups()
    else:
         d = dict(r.tags)
         m = d
         _id, _count, _len, _umi, _bc = r.qname, 'NA', len(r.query), d['XM'], d['XC']
    if m is not None:
        f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(_id,_umi,Seq(_umi).reverse_complement(),_bc,Seq(_bc).reverse_complement(),_len,_count))
    else: 
        m = rex_umi_only.match(r.description)
        _id, _count, _len, _umi = m.groups()
        _bc = 'NA'
        f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(_id,_umi,Seq(_umi).reverse_complement(),_bc,Seq(_bc).reverse_complement(),_len,_count))
f.close()
