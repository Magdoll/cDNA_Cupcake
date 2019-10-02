__author__ = 'etseng@pacb.com'

"""
To be used to evaluate phase switching by assuming the genome is already fully "phased".
For example, in Falcon-Unzip, which is pseudo-phased, there could be switching errors.
Since IsoPhase only uses the primary contig (or one of the haplotypes), a switching error shows up
 by seeing isoforms that have SNPs on the same allele being switched between different phases.

ex: in SNP 1, isoform 1 has GT=0|1
    in SNP 2, isoform 1 has GT=1|0
    in SNP 3, isoform 1 has GT=0|1

this would be a possible switching error. This is two switching events. SNP1->2 one event, 2->3 one event.
For now we will only consider the diploid case. Fail as soon as we see more than 2 alleles.
We will use phased.nopartial.cleaned.vcf as the input VCF.

Output is:

   <directory>, <chrom>, <strand>, <start>, <end>, <num_isoforms>, <num_switch_event>

We would get the loci (chr:start-end strand) from the config file in each directory.
"""


import os, sys, glob
from collections import defaultdict
import vcf


def read_config(config_filename):
    """
    pbid=PB.9923
    ref_chr=000104F_0
    ref_strand=-
    ref_start=3484015
    ref_end=3519963

    :return: (chr, start, end, strand)
    """
    for line in  open(config_filename):
        a,b = line.strip().split('=')
        if a == 'ref_chr': _chr = b
        elif a == 'ref_strand': _strand = b
        elif a == 'ref_start': _start = b
        elif a == 'ref_end': _end = b
    return _chr, _start, _end, _strand


def eval_isophase_phaseswitch(isophase_vcf, config_file, out_f, name='NA'):

    _chr, _start, _end, _strand = read_config(config_file)

    reader = vcf.VCFReader(open(isophase_vcf))
    # record the first SNP for each isoform
    prev = {} # sample -> CallData.GT (ex: '0|1')
    r = next(reader)
    for c in r.samples:
        prev[c.sample] = c.data.GT


    num_switch = 0

    for r in reader:
        for c in r.samples:
            if c.data.GT.find('|') == -1: continue # ignore those with just one allele
            a, b = c.data.GT.split('|')
            if a == b: continue # for now, ignore IsoPhase results that only uses one allele
            if prev[c.sample] != c.data.GT:
                num_switch += 1
            prev[c.sample] = c.data.GT

    out_f.write("{name}\t{chrom}\t{start}\t{end}\t{strand}\t{num_iso}\t{num_switch}\n".format(\
        name=name, chrom=_chr, start=_start, end=_end, strand=_strand,
        num_iso=len(r.samples), num_switch=num_switch))


def main_eval():

    out_f = open('evaled.isophase_phase_switches.txt', 'w')
    out_f.write("dir\tchrom\tstart\tend\tstrand\tnum_iso\tnum_phase_switch\n")
    dirs = glob.glob('by_loci/*size*/')

    for d1 in dirs:
        configfile = os.path.join(d1, 'config')
        vcffile = os.path.join(d1, 'phased.nopartial.cleaned.vcf')
        nosnp = os.path.join(d1, 'phased.nopartial.NO_SNPS_FOUND')
        nohap = os.path.join(d1, 'phased.nopartial.NO_HAPS_FOUND')

        if not os.path.exists(vcffile):
            # no SNP, just skip
            assert os.path.exists(nosnp) or os.path.exists(nohap)
            print("Skipping {0} because no SNPs found.".format(d1), file=sys.stderr)
        else:
            print("Evaluating {0}.".format(d1), file=sys.stderr)
            name = d1.split('/')[1]
            eval_isophase_phaseswitch(vcffile, configfile, out_f, name)
    out_f.close()


main_eval()