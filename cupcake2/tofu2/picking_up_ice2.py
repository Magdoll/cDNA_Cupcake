__author__ = 'etseng@pacb.com'

#!/usr/bin/env python
import os, sys
from cPickle import *
from pbcore.io.FastqIO import FastqReader, FastqWriter

import cupcake2.ice2.IceIterative2 as ice

from pbtranscript.ice.IceUtils import set_probqv_from_model

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

def ensure_pickle_goodness(pickle_filename, root_dir, fasta_files_to_add, fastq_files_to_add):
    """
    Old versions of IceIterative.write_pickle is missing some key/values.
    Add if needed.
    Return a good pickle object
    """
    a = load(open(pickle_filename))
    if a['fasta_filename'] != os.path.abspath(os.path.join(root_dir,'current.fasta')):
        raise Exception, "The pickle file {0} indicates that current.fasta is not being used. ICE2 likely did not finish to a point that could be picked up.".format(pickle_filename)

    a['newids'] = check_n_fix_newids(a)

    a['fasta_filenames_to_add'] = fasta_files_to_add.split(',')
    a['fastq_filenames_to_add'] = fastq_files_to_add.split(',')

    # add the paired fasta/fastq files
    if len(fasta_files_to_add) != len(fastq_files_to_add):
        raise Exception, "Number of fasta and fastq files to add is not equal!"

    # sanity check that the fasta/fastq files to add exist
    for file in fastq_files_to_add.split(','):
        if not os.path.exists(file):
            raise Exception, "{0} does not exist!".format(file)
    for file in fasta_files_to_add.split(','):
        if not os.path.exists(file):
            raise Exception, "{0} does not exist!".format(file)

    if 'root_dir' not in a:
        print >> sys.stderr, "Pickle {0} missing some key-values. Fixing it.".format(pickle_filename)
        a['root_dir'] = root_dir
        a['all_fasta_filename'] = a['all_fasta_fiilename']
        a['qv_prob_threshold'] = 0.03
        with open(pickle_filename + '.fixed', 'w') as f:
            dump(a, f)
        print >> sys.stderr, "Fixed pickle written to {0}.fixed".format(pickle_filename)
        return a, pickle_filename + '.fixed'
    else:
        # newid might have been fixed, STILL output pickle writing anyway
        with open(pickle_filename, 'w') as f:
            dump(a, f)
        return a, pickle_filename

def check_n_fix_newids(icec_obj):
    newids = icec_obj['newids']

    if len(newids) == 0:
        print >> sys.stderr, "newids is empty (probably a finalized run). set it."
        for k, v in icec_obj['d'].iteritems():
            if len(v) != 1:
                newids.add(k)
        print >> sys.stderr, "added {0} seqs to newids".format(len(newids))
    return newids


def pickup_icec_job(pickle_filename, root_dir, flnc_filename, fasta_files_to_add, fastq_files_to_add):
    icec_obj, icec_pickle_filename = ensure_pickle_goodness(pickle_filename, root_dir, fasta_files_to_add, fastq_files_to_add)
    probqv, msg = set_probqv_from_model()

    icec = ice.IceIterative2.from_pickle(icec_pickle_filename, probqv)

    # first must RE-RUN gcon to get all the proper refs
    icec.changes = set()
    icec.refs = {}
    icec.all_fasta_filename = flnc_filename
    todo = icec.uc.keys()
    print >> sys.stderr, "Re-run gcon for proper refs...."
    icec.run_gcon_parallel(todo)
    print >> sys.stderr, "Re-calculating cluster prob, just to be safe...."
    icec.calc_cluster_prob(True)
    print >> sys.stderr, "Sanity checking now...."
    icec.sanity_check_uc_refs()
    icec.ensure_probQV_newid_consistency()
    print >> sys.stderr, "Sanity check done. Resuming ICE job."
    icec.run()


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("pickle_filename", help="Last successful pickle (ex: clusterOut/output/input.split_00.fa.pickle)")
    parser.add_argument("--root_dir", default="clusterOut/", help="root dir (default: clusterOut/)")
    parser.add_argument("--flnc", default="isoseq_flnc.fasta", help="(default: isoseq_flnc.fasta)")
    parser.add_argument("--fasta_files_to_add", default=None, help="Comma-separated additional fasta files to add (default: None)")
    parser.add_argument("--fastq_files_to_add", default=None, help="Comma-separated additional fastq files to add (default: None)")

    args = parser.parse_args()

    pickup_icec_job(args.pickle_filename, args.root_dir, args.flnc, args.fasta_files_to_add, args.fasta_files_to_add)

