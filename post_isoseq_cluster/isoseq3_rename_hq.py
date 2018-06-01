#!/usr/bin/env python
import os, sys
from Bio import SeqIO

if not os.path.exists('unpolished.fasta'):
    print >> sys.stderr, "unpolished.fasta does not exist. Abort!"
    sys.exit(-1)
if not os.path.exists("polished.hq.fastq"):
    print >> sys.stderr, "polished.hq.fastq does not exist. Abort!"
    sys.exit(-1)

fl_count = {} # isoform index --> count
for r in SeqIO.parse(open('unpolished.fasta'),'fasta'):
    # id: isoform_0 full_length_coverage=2;length=12786
    i = r.id.split('_')[1]
    c = int(r.description.split('full_length_coverage=')[1].split(';')[0])
    fl_count[i] = c

h = open('hq.fastq', 'w')
reader = SeqIO.parse(open('polished.hq.fastq'),'fastq')
for r in reader:
    acc = 1 - sum(10**(-x/10.) for x in r.letter_annotations['phred_quality'][100:-30])/len(r.seq)
    i = r.id.split('/')[1]
    c = fl_count[i]
    r.id = "HQ_sampleX|cb1_c{0}/f{1}p0/{2}".format(i, c, len(r.seq))
    r.description = "{0} full_length_coverage={1};isoform_length={2};expected_accuracy={3:.4f}".format(r.id,c,len(r.seq),acc)
    SeqIO.write(r, h, 'fastq')
    
h.close()

print >> sys.stderr, "Output written to: hq.fastq"
