__author__ = 'etseng@pacb.com'

#!/usr/bin/env python
import os, sys
from pickle import *
from pbcore.io.FastqIO import FastqReader, FastqWriter

import cupcake2.ice2.IceIterative2 as ice

import pbtranscript.ice.ProbModel as pm

"""
Pickle requirements:

root_dir: root dir
changes: (ignore)
uc: (original uc from IceInit, ignore)
d: dict of seqid --> cid --> prob
refs: dict of cid --> gcon ref filename
newids: set of seqids being iterated (not frozen)
all_fasta_filename: isoseq_flnc.fasta
fasta_filename: (ignore, usually current.fasta)
qv_prob_threshold: usually set of 0.03
fasta_filenames_to_add: (not yet added files)

ice_opts: cupcake2.tofu2.ClusterOptions2.IceOptions2
sge_opts: cupcake2.tofu2.ClusterOptions2.SgeOptions2

"""

def ensure_pickle_goodness(pickle_filename, root_dir, flnc_filename, fasta_files_to_add, fastq_files_to_add):
    """
    Old versions of IceIterative.write_pickle is missing some key/values.
    Add if needed.
    Return a good pickle object
    """
    a = load(open(pickle_filename))
    if a['fasta_filename'] != os.path.abspath(os.path.join(root_dir,'current.fasta')):
        raise Exception("The pickle file {0} indicates that current.fasta is not being used. ICE2 likely did not finish to a point that could be picked up.".format(pickle_filename))

    a['newids'] = check_n_fix_newids(a)

    a['fasta_filenames_to_add'] = fasta_files_to_add.split(',')
    a['fastq_filenames_to_add'] = fastq_files_to_add.split(',')
    a['all_fasta_filename'] = flnc_filename

    # add the paired fasta/fastq files
    if len(fasta_files_to_add) != len(fastq_files_to_add):
        raise Exception("Number of fasta and fastq files to add is not equal!")

    # sanity check that the fasta/fastq files to add exist
    for file in fastq_files_to_add.split(','):
        if not os.path.exists(file):
            raise Exception("{0} does not exist!".format(file))
    for file in fasta_files_to_add.split(','):
        if not os.path.exists(file):
            raise Exception("{0} does not exist!".format(file))

    if 'root_dir' not in a:
        print("Pickle {0} missing some key-values. Fixing it.".format(pickle_filename), file=sys.stderr)
        a['root_dir'] = root_dir
        a['qv_prob_threshold'] = 0.03
        with open(pickle_filename + '.fixed', 'w') as f:
            dump(a, f)
        print("Fixed pickle written to {0}.fixed".format(pickle_filename), file=sys.stderr)
        return a, pickle_filename + '.fixed'
    else:
        # newid might have been fixed, STILL output pickle writing anyway
        with open(pickle_filename, 'w') as f:
            dump(a, f)
        return a, pickle_filename

def check_n_fix_newids(icec_obj):
    newids = icec_obj['newids']

    if len(newids) == 0:
        print("newids is empty (probably a finalized run). set it.", file=sys.stderr)
        for k, v in icec_obj['d'].items():
            if len(v) != 1:
                newids.add(k)
        print("added {0} seqs to newids".format(len(newids)), file=sys.stderr)
    return newids


def make_current_fastq(icec_obj, flnc_filename, root_dir):
    """
    current fasta will consists of all ids

    however --- if this was a already finished run and we are adding more input,
        then newids is empty, in this case we set newids = everything that
        has no affiliation or more than one affiliated cluster in d
    """
    with FastqWriter(os.path.join(root_dir, 'current.fastq')) as f:
        for r in FastqReader(flnc_filename):
            f.writeRecord(r)

def pickup_icec_job(pickle_filename, root_dir, flnc_filename, flnc_fq, fasta_files_to_add, fastq_files_to_add):
    icec_obj, icec_pickle_filename = ensure_pickle_goodness(pickle_filename, root_dir, flnc_filename, fasta_files_to_add, fastq_files_to_add)
    make_current_fastq(icec_obj, flnc_fq, root_dir)
    print("Reading QV information....", file=sys.stderr)
    probqv = pm.ProbFromFastq(os.path.join(root_dir,'current.fastq'))

    icec = ice.IceIterative2.from_pickle(icec_pickle_filename, probqv)

    # first must RE-RUN gcon to get all the proper refs
    icec.changes = set()
    icec.refs = {}
    icec.all_fasta_filename = flnc_filename
    icec.fasta_filename = 'current.fasta'
    icec.fastq_filename = 'current.fastq'
    icec.fasta_filenames_to_add = fasta_files_to_add.split(',')
    icec.fastq_filenames_to_add = fastq_files_to_add.split(',')

    todo = list(icec.uc.keys())
    print("Re-run gcon for proper refs....", file=sys.stderr)
    icec.run_gcon_parallel(todo)
    print("Re-calculating cluster prob, just to be safe....", file=sys.stderr)
    icec.calc_cluster_prob(True)
    print("Sanity checking now....", file=sys.stderr)
    icec.sanity_check_uc_refs()
    icec.ensure_probQV_newid_consistency()
    print("Sanity check done. Resuming ICE job.", file=sys.stderr)
    icec.run()


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("pickle_filename", help="Last successful pickle (ex: clusterOut/output/input.split_00.fa.pickle)")
    parser.add_argument("--root_dir", default="clusterOut/", help="root dir (default: clusterOut/)")
    parser.add_argument("--flnc_fa", default="isoseq_flnc.fasta", help="(default: isoseq_flnc.fasta)")
    parser.add_argument("--flnc_fq", default="isoseq_flnc.fastq", help="(default: isoseq_flnc.fastq)")
    parser.add_argument("--fasta_files_to_add", default=None, help="Comma-separated additional fasta files to add (default: None)")
    parser.add_argument("--fastq_files_to_add", default=None, help="Comma-separated additional fastq files to add (default: None)")

    args = parser.parse_args()

    pickup_icec_job(args.pickle_filename, args.root_dir, args.flnc_fa, args.flnc_fq, args.fasta_files_to_add, args.fastq_files_to_add)

