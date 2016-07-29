#!/usr/bin/env python
import os, sys
from collections import defaultdict
from Bio import SeqIO


fasta_filename = sys.argv[1]

d = defaultdict(lambda: 0)
for r in SeqIO.parse(open(fasta_filename),'fasta'):
    d[r.id.split('/')[0]] += 1

for k,v in d.iteritems(): print k,v
