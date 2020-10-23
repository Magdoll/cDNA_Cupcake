#!/usr/bin/env python
__author__ = 'etseng@pacb.com'

"""
Similar to chain_samples.py, but specifically for chaining fusions.

Difference with regular chaining:

1. A fusion consists of several transcript records, denoted by PB.X.1, PB.X.2, ... where PB.X is one fusion
   (must always consider possibility of a fusion consisting of 2+ loci)

2. For two fusions from two samples to be "equivalent", every loci (.1, .2, ..) must agree

"""

from cupcake.tofu.counting.chain_samples import read_config, read_count_info


import os, sys, glob, shutil
from csv import DictReader
from collections import defaultdict
from Bio import SeqIO
import cupcake.io.GFF as GFF
from cupcake.tofu.counting import combine_abundance_across_samples as sp

def sample_sanity_check(group_filename, gff_filename, count_filename, fastq_filename=None):
    """
    Double check that the formats are expected and all PBIDs are concordant across the files
    :return: raise Exception if sanity check failed
    """
    print("Sanity checking. Retrieving PBIDs from {0},{1},{2}...".format(\
        group_filename, gff_filename, count_filename), file=sys.stderr)
    ids1 = [line.strip().split()[0] for line in open(group_filename)]
    ids2 = [fusion_id for fusion_id,rs in GFF.collapseGFFFusionReader(gff_filename)]
    f = open(count_filename)
    while True:
        # advance through the headers which start with #
        cur = f.tell()
        if not f.readline().startswith('#') or f.tell() == cur:  # first non-# seen or EOF
            f.seek(cur)
            break
    ids3 = [r['pbid'] for r in DictReader(f, delimiter='\t')]
    if len(set(ids2).difference(ids1))>0 or len(set(ids2).difference(ids3))>0:
        raise Exception("Sanity check failed! Please make sure the PBIDs listed in {1} are also in {0} and {2}".format(\
            group_filename, gff_filename, count_filename))

    if fastq_filename is not None:
        ids4 = [r.id.split('|')[0] for r in SeqIO.parse(open(fastq_filename), 'fastq')]
        if len(set(ids2).difference(ids4))>0:
            raise Exception("Sanity check failed! Please make sure the PBIDs listed in {1} are also in {0}".format(\
                fastq_filename, gff_filename))


def chain_fusion_samples(dirs, names, group_filename, gff_filename, count_filename, field_to_use='count_fl', fuzzy_junction=0, fastq_filename=None):

    for d in dirs.values():
        sample_sanity_check(os.path.join(d, group_filename),\
                            os.path.join(d, gff_filename),\
                            os.path.join(d, count_filename),\
                            os.path.join(d, fastq_filename) if fastq_filename is not None else None)

    count_header, count_info = read_count_info(count_filename, dirs, field_to_use)

    # some names may already start with "tmp_" which means they are intermediate results that have already been chained
    # find the first non "tmp_" and start from there
    if names[0].startswith('tmp_'):
        chain = []
        for start_i,name in enumerate(names):
            if name.startswith('tmp_'):
                chain.append(name[4:])
            else:
                break
        # start_i, name now points at the first "non-tmp" sample
        # we want to go to the last tmp_ sample and read it
        name = names[start_i-1][4:] # this is the last tmp_ sample, let's read it
        o = sp.MegaPBTreeFusion('tmp_'+name+'.gff', 'tmp_'+name+'.group.txt', self_prefix='tmp_'+name, \
                        internal_fuzzy_max_dist=fuzzy_junction, \
                        fastq_filename='tmp_'+name+'.rep.fq' if fastq_filename is not None else None)
        #chain.append(name) # no need, already done above
    else: # everything is new, start fresh
        name = names[0]
        d = dirs[name]
        chain = [name]
        o = sp.MegaPBTreeFusion(os.path.join(d, gff_filename), os.path.join(d, group_filename), \
                        self_prefix=name, internal_fuzzy_max_dist=fuzzy_junction, \
                        fastq_filename=os.path.join(d, fastq_filename) if fastq_filename is not None else None)
        start_i = 1

    for name in names[start_i:]:
        assert not name.startswith('tmp_')
        d = dirs[name]
        o.add_sample(os.path.join(d, gff_filename), os.path.join(d, group_filename), \
                     sample_prefix=name, output_prefix='tmp_'+name, \
                     fastq_filename=os.path.join(d, fastq_filename) if fastq_filename is not None else None)
        o = sp.MegaPBTreeFusion('tmp_'+name+'.gff', 'tmp_'+name+'.group.txt', self_prefix='tmp_'+name, \
                          internal_fuzzy_max_dist=fuzzy_junction, \
                          fastq_filename='tmp_'+name+'.rep.fq' if fastq_filename is not None else None)
        chain.append(name)

    # now recursively chain back by looking at mega_info.txt!!!
    d = {} # ex: (tmp_1009, PB.1.1) --> mega info dict
    for c in chain[1:]:
        for r in DictReader(open('tmp_' + c + '.mega_info.txt'),delimiter='\t'):
            d['tmp_'+c, r['pbid']] = r

    f1 = open('all_samples.chained_ids.txt', 'w')
    f1.write("superPBID")
    f2 = open('all_samples.chained_count.txt', 'w')
    f2.write("superPBID")
    for c in chain:
        f1.write('\t' + c)
        f2.write('\t' + c)
    f1.write('\n')
    f2.write('\n')

    reader = DictReader(open('tmp_' + chain[-1] + '.mega_info.txt'),delimiter='\t')
    for r in reader:
        saw_NA = False
        r0 = r
        answer = defaultdict(lambda: 'NA') # ex: 1009 --> PB.1.1
        answer2 = defaultdict(lambda: 'NA') # ex: 1009 --> count
        answer[chain[-1]] = r[chain[-1]]
        if r[chain[-1]] !='NA':
            answer2[chain[-1]] = count_info[chain[-1], answer[chain[-1]]]
        for c in chain[::-1][1:-1]:  # the first sample does not have tmp_, because it's not a chain
            if r['tmp_'+c] == 'NA':
                saw_NA = True
                break
            else:
                r2 = d['tmp_'+c, r['tmp_'+c]]
                answer[c] = r2[c]
                if answer[c] != 'NA':
                    answer2[c] = count_info[c, answer[c]]
                r = r2
        if not saw_NA:
            answer[chain[0]] = r[chain[0]]
            if answer[chain[0]] !='NA':
                answer2[chain[0]] = count_info[chain[0], answer[chain[0]]]
        f1.write(r0['pbid'])
        f2.write(r0['pbid'])
        for c in chain:
            f1.write("\t" + answer[c]) # each tissue still share the same PB id
            f2.write("\t" + str(answer2[c]))
        f1.write('\n')
        f2.write('\n')
    f1.close()
    f2.close()

    shutil.copyfile('tmp_' + chain[-1] + '.gff', 'all_samples.chained.gff')
    if fastq_filename is not None:
        shutil.copyfile('tmp_' + chain[-1] + '.rep.fq', 'all_samples.chained.rep.fq')

    print("Chained output written to:", file=sys.stderr)
    print("all_samples.chained.gff", file=sys.stderr)
    print(f1.name, file=sys.stderr)
    print(f2.name, file=sys.stderr)
    if fastq_filename is not None:
        print("all_samples.chained.rep.fq", file=sys.stderr)


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("config_file")
    parser.add_argument("field_to_use", choices=['norm_fl', 'count_fl'], default='count_fl', help="Which count field to use for chained sample (default: count_fl)")
    parser.add_argument("--fuzzy_junction", default=5, type=int, help="Max allowed distance in junction to be considered identical (default: 5 bp)")
    #parser.add_argument("--allow_5merge", action="store_true", default=False, help="Allow 5' truncated transcripts (default: off)" )
    args = parser.parse_args()

    sample_dirs, sample_names, group_filename, gff_filename, count_filename, fastq_filename = read_config(args.config_file)
    chain_fusion_samples(sample_dirs, sample_names, group_filename, gff_filename, count_filename, args.field_to_use, args.fuzzy_junction, fastq_filename)


