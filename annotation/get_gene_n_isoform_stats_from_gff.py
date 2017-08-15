#!/usr/bin/env python
import cupcake.io.GFF as GFF
from collections import defaultdict
import numpy as np

d = defaultdict(lambda: [])
for r in GFF.collapseGFFReader('hq_isoforms.fastq.no5merge.collapsed.filtered.gff'):
    d[r.seqid.split('.')[1]].append(r.seqid)
    
p = np.array([len(v) for v in d.itervalues()])

print "Number of loci:", len(d)
print "Number of isoforms:", sum(p)
print "Avg. number of isoforms per loci:", np.mean(p)
f = open('hq_isoforms.fastq.no5merge.collapsed.filtered.isoform_per_loci.txt', 'w')
f.write("loci\tnum_isoform\n")
for k,v in d.iteritems(): f.write("PB.{0}\t{1}\n".format(k, len(v)))
f.close()
