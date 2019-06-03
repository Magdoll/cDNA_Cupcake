import os
import vcf
import glob
from collections import defaultdict

dirs = glob.glob('by_loci/*size*')
snps_by_chrom = defaultdict(lambda: [])

for d in dirs:
    filename = os.path.join(d, 'phased.partial.vcf')
    if not os.path.exists(filename): continue
    reader = vcf.VCFReader(open(filename))
    reader.samples = []
    for r in reader: snps_by_chrom[r.CHROM].append((r.POS, r))

keys = snps_by_chrom.keys()
keys.sort()    
        
f = vcf.Writer(open('IsoSeq_IsoPhase.vcf', 'w'), reader)
for k in keys:
    v = snps_by_chrom[k]
    v.sort(key=lambda x: x[0])
    for pos,rec in v:
        rec.samples = []
        f.write_record(rec)
        
f.close()
