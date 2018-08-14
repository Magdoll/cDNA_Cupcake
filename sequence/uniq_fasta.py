#!/usr/bin/env python
import os, sys
from Bio import SeqIO

input = sys.argv[1]
output = input[:input.rfind('.')]+'.dedup.fasta'

seen = set()
f = open(output, 'w')
for r in SeqIO.parse(open(input),'fasta'):
    if r.id not in seen: f.write(">{0}\n{1}\n".format(r.description, r.seq))
    seen.add(r.id)

print >> sys.stderr, "De-dup fasta written to: {0}".format(output)
