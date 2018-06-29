# uncompyle6 version 3.2.1
# Python bytecode 2.7 (62211)
# Decompiled from: Python 2.7.14 |Anaconda custom (64-bit)| (default, Dec  7 2017, 17:05:42) 
# [GCC 7.2.0]
# Embedded file name: evaluate_snp_with_genome.py
# Compiled at: 2018-01-02 12:49:59
__author__ = 'etseng@pacb.com'
import os, sys, glob
from collections import defaultdict
import vcf
import phasing.io.SAMMPileUpReader as sp

def read_fake_mapping(fake_genome_mapping_filename):
    """
    :param fake_genome_mapping_filename: fake genome mapping file
    :return: dict of fake position -> (chr, ref position)
    """
    fake_map = {}
    with open(fake_genome_mapping_filename) as (f):
        for line in f:
            fake_pos, ref_chr, ref_pos = line.strip().split(',')
            fake_map[int(fake_pos)] = (ref_chr, int(ref_pos))

    return fake_map


def get_positions_to_recover(fake_genome_mapping_filename, mpileup_filename, genome_snp, min_cov=40):
    """
    Return a list of genomic positions that isophase should call SNPs for. Must be:
    
    1. covered (based on ccs.mpileup) with at least MIN_COVERAGE (default: 40)
    2. is a SNP (based on genome VCF)
    """
    fake_map = read_fake_mapping(fake_genome_mapping_filename)
    good_positions = []
    for r in sp.MPileUpReader(mpileup_filename):
        if r.nCov >= min_cov:
            genome_chr, genome_pos = fake_map[r.pos]
            genome_chr = genome_chr.split('|')[0]
            if genome_chr in genome_snp and genome_pos in genome_snp[genome_chr] and genome_snp[genome_chr][genome_pos].is_snp:
                good_positions.append((genome_chr, genome_pos))

    return good_positions


def eval_isophase(isophase_vcf, genome_snp, good_positions, out_f, name='NA'):
    for r in vcf.VCFReader(open(isophase_vcf)):
        r.CHROM = r.CHROM.split('|')[0]
        if (r.CHROM, r.POS) not in good_positions:
            alt_g = 'NA'
            in_g = 'N'
            in_i = 'Y'
        else:
            alt_g = genome_snp[r.CHROM][r.POS].ALT[0]
            in_g = 'Y'
            in_i = 'Y'
            good_positions.remove((r.CHROM, r.POS))
        out_f.write(('{name}\t{chrom}\t{pos}\t{ref}\t{alt_g}\t{alt_i}\t{in_g}\t{in_i}\n').format(name=name, chrom=r.CHROM, pos=r.POS, ref=r.REF, alt_g=alt_g, alt_i=r.ALT[0], in_g=in_g, in_i=in_i))

    for chrom, pos in good_positions:
        r = genome_snp[chrom][pos]
        out_f.write(('{name}\t{chrom}\t{pos}\t{ref}\t{alt_g}\tNA\tY\tN\n').format(name=name, chrom=chrom, pos=pos, ref=r.REF, alt_g=r.ALT[0]))


def main_brangus(unzip_snps=None):
    if unzip_snps is None:
        unzip_snps = defaultdict(lambda : {})
        for r in vcf.VCFReader(open('Brangus.unzip.vcf')):
            unzip_snps[r.CHROM][r.POS] = r

    print >> sys.stderr, 'Finished reading Brangus.unzip.vcf.'
    out_f = open('evaled.isophase.txt', 'w')
    out_f.write('dir\tchrom\tpos\tref\talt_g\talt_i\tin_g\tin_i\n')
    dirs = glob.glob('by_loci/*size*/')
    for d1 in dirs:
        mpileup = os.path.join(d1, 'ccs.mpileup')
        mapfile = os.path.join(d1, 'fake.mapping.txt')
        vcffile = os.path.join(d1, 'phased.partial.vcf')
        nosnp = os.path.join(d1, 'phased.partial.NO_SNPS_FOUND')
        if not os.path.exists(vcffile):
            assert os.path.exists(nosnp)
            print >> sys.stderr, ('Skipping {0} because no SNPs found.').format(d1)
        else:
            print >> sys.stderr, ('Evaluating {0}.').format(d1)
            good_positions = get_positions_to_recover(mapfile, mpileup, unzip_snps, min_cov=40)
            name = d1.split('/')[1]
            eval_isophase(vcffile, unzip_snps, good_positions, out_f, name)

    out_f.close()
    return

def main_maize(ki11_snps=None, dirs=None):
    if ki11_snps is None:
        ki11_snps = defaultdict(lambda : {}) # chrom -> pos -> VCF record
        for r in vcf.VCFReader(open('B73Ki11.q20.vcf')):
            ki11_snps[r.CHROM][r.POS] = r

    print >> sys.stderr, 'Finished reading B73Ki11.q20.vcf.'
    out_f = open('evaled.isophase.txt', 'w')
    out_f.write('dir\tchrom\tpos\tref\talt_g\talt_i\tin_g\tin_i\n')
    if dirs is None: dirs = glob.glob('by_loci/*size*/')
    for d1 in dirs:
        mpileup = os.path.join(d1, 'ccs.mpileup')
        mapfile = os.path.join(d1, 'fake.mapping.txt')
        vcffile = os.path.join(d1, 'phased.partial.vcf')
        nosnp = os.path.join(d1, 'phased.partial.NO_SNPS_FOUND')
        if not os.path.exists(vcffile):
            assert os.path.exists(nosnp)
            print >> sys.stderr, ('Skipping {0} because no SNPs found.').format(d1)
        else:
            print >> sys.stderr, ('Evaluating {0}.').format(d1)
            good_positions = get_positions_to_recover(mapfile, mpileup, ki11_snps, min_cov=30) # use lower min cov here becuz a few close cases where BQ filtering lowered cov
            name = d1.split('/')[1]
            eval_isophase(vcffile, ki11_snps, good_positions, out_f, name)

    out_f.close()
    return ki11_snps

if __name__ == "__main__":
    from csv import DictReader
    dirs = ['by_loci/'+r['locus'] for r in DictReader(open('evaled_isophase.demux_hap_count.txt'),delimiter='\t')]
    main_maize(None, dirs)
