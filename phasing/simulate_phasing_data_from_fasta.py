#!/usr/bin/env python
__author__ = 'etseng@pacb.com'

"""
Input: fasta, each sequence is an individual gene

given an input transcript sequence ---
-> randomly select N positions to have a SNP, N = <length> / 300
-> generate two haplotypes, one is the original, one is mutated
(store this in Haplotypes)
-> for each haplotype, simulate 100X sequences at 0%, 1%, 2%, 3% (each is a diff output fasta)

"""
import os, sys, random
from collections import Counter
from Bio import SeqIO
from simulate import sim_seq
from phasing.io.VariantPhaser import Haplotypes

SNP_INTERVAL = 300  # 1 snp per 300 bp
base_choices = {'A': ('T','C','G'), 'T': ('A','C','G'), 'C': ('A','T','G'), 'G':('A','T','C')}

def simulate_phasing_data(seq0, err_sub, ploidity, copies, write_fastq=False):
    """
    :param seq0: transcript sequence
    :param err_sub: error prob
    :param ploidity: ploidity (how many alleles)
    :param copies: list of how many copies to simulate per allele   ex: [10, 10] means 10 for each allele
    """
    n = len(seq0)
    var_positions = random.sample(list(range(n)), n/SNP_INTERVAL)
    var_positions.sort()

    new_seqs = [seq0]
    for ignore in range(ploidity-1):
        num_tries = 0
        while num_tries < 10:  # after 10 attemps, give up simulating enough divergent seqs (this happens for very short haps)
            num_tries += 1
            new_seq = list(seq0)
            for p in var_positions:
                new_seq[p] = random.choice(base_choices[seq0[p]])
            new_seq = "".join(new_seq)
            if new_seq not in new_seqs:
                new_seqs.append("".join(new_seq))
                break


    count_of_vars_by_pos = {}
    for p in var_positions:
        count_of_vars_by_pos[p] = Counter() # base --> count of base at pos
    for i,seq in enumerate(new_seqs):
        for p in var_positions:
            count_of_vars_by_pos[p][seq[p]] += copies[i]

    # 1. write out fake.fasta
    f = open('fake.fasta', 'w')
    f.write(">fake\n{0}\n".format(seq0))
    f.close()
    # 2. write fake.mapping.txt
    f = open('fake.mapping.txt', 'w')
    for i in range(n): f.write("{0},fake,{0}\n".format(i))
    f.close()

    # simulate CCS reads
    if write_fastq:
        f = open('ccs.fastq', 'w')
    else:
        f = open('ccs.fasta', 'w')
    f2 = open('fake.read_stat.txt', 'w')
    f2.write("id\tlength\tis_fl\tstat\tpbid\n")
    profile = [err_sub, err_sub, err_sub, 1.]
    for i,copy in enumerate(copies):
        if i >= len(new_seqs): break
        for c in range(copy):
            cur_seq, cur_qv = sim_seq(new_seqs[i], profile)
            cur_id = "hap{0}_{1}".format(i+1, c)
            if write_fastq:
                f.write("@{0}\n{1}\n".format(cur_id, cur_seq))
                f.write("+\n{0}\n".format(cur_qv))
            else:
                f.write(">{0}\n{1}\n".format(cur_id, cur_seq))
            f2.write("{0}\t{1}\tY\tunique\tPB.0.0\n".format(cur_id, n))

    f.close()
    f2.close()

    ref_at_pos = dict((p,seq0[p]) for p in var_positions)
    hap_obj = Haplotypes(var_positions, ref_at_pos, count_of_vars_by_pos)
    for seq in new_seqs:
        hap = "".join(seq[p] for p in var_positions)
        hap_obj.match_or_add_haplotype(hap)

    hap_obj.get_haplotype_vcf_assignment()
    hap_obj.write_haplotype_to_vcf('fake.mapping.txt', {}, 'fake.answer')

if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("fasta_filename", help="Fasta file from which to simulate phasing data from.")
    parser.add_argument("-p", "--ploidity", type=int, default=2)
    parser.add_argument("--err_sub", type=float, required=True)
    parser.add_argument("--copies", required=True)
    parser.add_argument("--write_fastq", action="store_true", default=False)


    args = parser.parse_args()

    assert 2 <= args.ploidity <= 6

    copies = list(map(int, args.copies.split(',')))
    assert len(copies) == args.ploidity

    for r in SeqIO.parse(open(args.fasta_filename), 'fasta'):
        d2 = r.id.split('|')[0]
        print("making {0}".format(d2), file=sys.stderr)
        os.makedirs(d2)
        os.chdir(d2)
        simulate_phasing_data(r.seq.tostring(), args.err_sub, args.ploidity, copies, args.write_fastq)
        os.chdir('../')
