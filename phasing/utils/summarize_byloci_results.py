import os, sys, glob
import vcf
from csv import DictWriter, DictReader
from Bio import SeqIO

FIELDS = ['locus', 'size', 'num_snp', 'num_hap_nopartial', 'num_hap_withpartial']

dirs = glob.glob('by_loci/*size*')

f = open('summarized.isophase_results.txt', 'w')
writer = DictWriter(f, FIELDS, delimiter='\t')
writer.writeheader()

for d in dirs:
    size = 0
    for r in SeqIO.parse(open(os.path.join(d, "ccs.fasta")), 'fasta'): size += 1

    rec = {'locus': d, 'size': size}

    if os.path.exists(os.path.join(d, 'phased.nopartial.NO_SNPS_FOUND')):
        rec['num_snp'] = 0
        rec['num_hap_nopartial'] = 0
        rec['num_hap_withpartial'] = 0
    elif not os.path.exists(os.path.join(d, 'phased.partial.vcf')):
        print("WARNING: cannot find {0}/phased.partial.vcf. IGNORE".format(d))
    else:
        rec['num_snp'] = len([x for x in vcf.VCFReader(open(os.path.join(d, 'phased.partial.vcf')))])
        if os.path.exists(os.path.join(d, 'phased.nopartial.NO_HAPS_FOUND')):
            rec['num_hap_nopartial'] = 0
            rec['num_hap_withpartial'] = 0
        else:
            file1 = os.path.join(d, 'phased.nopartial.cleaned.human_readable.txt')
            file2 = os.path.join(d, 'phased.partial.cleaned.human_readable.txt')
            h1 = open(file1)
            h1.readline() # skip header
            h2 = open(file2)
            h2.readline() # skip header
            rec['num_hap_nopartial'] = len([r for r in DictReader(h1, delimiter='\t')])
            rec['num_hap_withpartial'] = len([r for r in DictReader(h2, delimiter='\t')])
    writer.writerow(rec)
f.close()

print("Summarized results of by_loci/<dirs> to {0}.".format(f.name), file=sys.stderr)
