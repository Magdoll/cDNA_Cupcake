__author__ = 'etseng@pacb.com'

#!/usr/bin/env python
import os, sys, glob, shutil
from csv import DictReader
from collections import defaultdict
from Bio import SeqIO
from cupcake.io import GFF
from cupcake.tofu.counting import combine_abundance_across_samples as sp

def sample_sanity_check(group_filename, gff_filename, count_filename, fastq_filename=None):
    """
    Double check that the formats are expected and all PBIDs are concordant across the files
    :return: raise Exception if sanity check failed
    """
    print >> sys.stderr, "Sanity checking. Retrieving PBIDs from {0},{1},{2}...".format(\
        group_filename, gff_filename, count_filename)
    ids1 = [line.strip().split()[0] for line in open(group_filename)]
    ids2 = [r.seqid for r in GFF.collapseGFFReader(gff_filename)]
    f = open(count_filename)
    for i in xrange(14): f.readline() # just through the header
    ids3 = [r['pbid'] for r in DictReader(f, delimiter='\t')]
    if set(ids1)!=set(ids2) or set(ids2)!=set(ids3):
        raise Exception, "Sanity check failed! Please make sure the PBIDs are consistent between: {0],{1},{2}".format(\
            group_filename, gff_filename, count_filename)

    if fastq_filename is not None:
        ids4 = [r.id.split('|')[0] for r in SeqIO.parse(open(fastq_filename), 'fastq')]
        if set(ids1)!=set(ids4):
            raise Exception, "Sanity check failed! {0} IDs don't match with the group file {1}".format(\
                fastq_filename, group_filename)


def read_config(filename):
    """
    SAMPLE=<name>;<path>

    must also have
    GROUP_FILENAME=
    GFF_FILENAME=
    COUNT_FILENAME=

    optional:
    FASTQ_FILENAME=
    """
    sample_dirs = {}
    sample_names = []
    group_filename, gff_filename, count_filename = None, None, None
    fastq_filename = None

    with open(filename) as f:
        for line in f:
            if line.startswith('SAMPLE='):
                name, path = line.strip()[7:].split(';')
                sample_dirs[name] = os.path.abspath(path)
                sample_names.append(name)
            elif line.startswith('GROUP_FILENAME='):
                group_filename = line.strip()[len('GROUP_FILENAME='):]
            elif line.startswith('GFF_FILENAME='):
                gff_filename = line.strip()[len('GFF_FILENAME='):]
            elif line.startswith('COUNT_FILENAME='):
                count_filename = line.strip()[len('COUNT_FILENAME='):]
            elif line.startswith('FASTQ_FILENAME='):
                fastq_filename = line.strip()[len('FASTQ_FILENAME='):]

    if group_filename is None:
        raise Exception, "Expected GROUP_FILENAME= but not in config file {0}! Abort.".format(filename)
    if count_filename is None:
        raise Exception, "Expected COUNT_FILENAME= but not in config file {0}! Abort.".format(filename)
    if gff_filename is None:
        raise Exception, "Expected GFF_FILENAME= but not in config file {0}! Abort.".format(filename)

    if len(sample_names) == 0:
        print >> sys.stderr, "No samples given. Exit."
        sys.exit(-1)

    return sample_dirs, sample_names, group_filename, gff_filename, count_filename, fastq_filename


def chain_samples(dirs, names, group_filename, gff_filename, count_filename, field_to_use='norm_nfl', fuzzy_junction=0, allow_5merge=False, fastq_filename=None):

    for d in dirs:
        sample_sanity_check(os.path.join(d, group_filename),\
                            os.path.join(d, gff_filename),\
                            os.path.join(d, count_filename),\
                            os.path.join(d, fastq_filename) if fastq_filename is not None else None)

    count_info = {} # key: (sample, PB.1.1) --> count
    for name, d in dirs.iteritems():
        f = open(os.path.join(d, count_filename))
        while True:
            cur = f.tell()
            if not f.readline().startswith('#'): break
        f.seek(cur)
        for r in DictReader(f, delimiter='\t'):
            count_info[name, r['pbid']] = r[field_to_use]

    name = names[0]
    d = dirs[name]
    chain = [name]

    o = sp.MegaPBTree(os.path.join(d, gff_filename), os.path.join(d, group_filename), \
                      self_prefix=name, internal_fuzzy_max_dist=fuzzy_junction, \
                      allow_5merge=allow_5merge, \
                      fastq_filename=os.path.join(d, fastq_filename))
    for name in names[1:]:
        d = dirs[name]
        o.add_sample(os.path.join(d, gff_filename), os.path.join(d, group_filename), \
                     sample_prefix=name, output_prefix='tmp_'+name, \
                     fastq_filename=os.path.join(d, fastq_filename))
        o = sp.MegaPBTree('tmp_'+name+'.gff', 'tmp_'+name+'.group.txt', self_prefix='tmp_'+name, \
                          internal_fuzzy_max_dist=fuzzy_junction, \
                          allow_5merge=allow_5merge, \
                          fastq_filename='tmp_'+name+'.rep.fq')
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

    print >> sys.stderr, "Chained output written to:"
    print >> sys.stderr, "all_samples.chained.gff"
    print >> sys.stderr, f1.name
    print >> sys.stderr, f2.name
    if fastq_filename is not None:
        print >> sys.stderr, "all_samples.chained.rep.fq"


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("config_file")
    parser.add_argument("field_to_use", choices=['norm_fl', 'norm_nfl', 'norm_nfl_amb', 'count_fl', 'count_nfl', 'count_nfl_amb'], default='norm_nfl', help="Which count field to use for chained sample (default: norm_nfl)")
    parser.add_argument("--fuzzy_junction", default=5, type=int, help="Max allowed distance in junction to be considered identical (default: 5 bp)")
    parser.add_argument("--allow_5merge", action="store_true", default=False, help="Allow 5' truncated transcripts (default: off)" )
    args = parser.parse_args()

    sample_dirs, sample_names, group_filename, gff_filename, count_filename, fastq_filename = read_config(args.config_file)
    chain_samples(sample_dirs, sample_names, group_filename, gff_filename, count_filename, args.field_to_use, args.fuzzy_junction, args.allow_5merge, fastq_filename)

