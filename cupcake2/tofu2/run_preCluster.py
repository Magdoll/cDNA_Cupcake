

import os, sys, pdb, subprocess
from cPickle import dump
from Bio import SeqIO
from cupcake2.ice2.preCluster import preClusterSet
from cupcake2.io.minimapIO import MiniReader
from cupcake.io.SeqReaders import LazyFastaReader
import cupcake2.ice2.preClusterProcess as sp
import cupcake2.ice2.AlignerRunners as ar
import cupcake2.io.FileIO as FileIO
from collections import defaultdict


def detect_PCR_chimeras(orphans, fasta_d, window_size=20, min_T_count=16):
    """
    Remove from orphans, all sequences that have stretch of polyTs in the front
    Return: candidates that fit the PCR chimera criteria
    """
    result = []
    for _id in orphans:
        r = fasta_d[_id]
        if r.seq[:window_size].count('T') >= min_T_count:
            result.append(r.id)
    return result


def add_batch(batch_index, pCS, orphans, fasta_d):
    """
    1. align batch<i>.fasta against seed<i>.S.fasta, process -> write remains to batch<i>.remains.fasta
    2. align batch<i>.remains.fasta against seed<i>.orphans.fasta -> write remains to batch<i>.remains2.fasta
    3. self align batch<i>.remains2.fasta -> combine remains+orphans to new orphans
    4. write out seed<i+1>.S.fasta and seed<i+1>.orphans.fasta

    """
    cur_file = "batch{0}.fasta".format(batch_index)
    seqids = set([r.id for r in SeqIO.parse(open(cur_file), 'fasta')])
    o = ar.run_minimap(cur_file, "seed{0}.S.fasta".format(batch_index))
    print >> sys.stderr, "processing", o
    pCS, remains = sp.process_align_to_pCS(o, seqids, pCS, MiniReader)
    print >> sys.stderr, "pCS: {0}, tucked: {1}, orphans: {2}, remains: {3}".format( \
        len(pCS.S), sum(v == 'T' for v in pCS.seq_stat.itervalues()), len(orphans), len(remains))
    # write batch<i>.remains.fasta
    cur_file = "batch{0}.remains.fasta".format(batch_index)
    FileIO.write_seqids_to_fasta(remains, cur_file, fasta_d)
    o = ar.run_minimap(cur_file, "seed{0}.orphans.fasta".format(batch_index))
    print >> sys.stderr, "processing", o
    pCS, orphans, remains = sp.process_align_to_orphan(o, remains, orphans, pCS, MiniReader)
    print >> sys.stderr, "pCS: {0}, tucked: {1}, orphans: {2}, remains: {3}".format( \
        len(pCS.S), sum(v == 'T' for v in pCS.seq_stat.itervalues()), len(orphans), len(remains))
    # write batch<i>.remains2.fasta and self align
    cur_file = "batch{0}.remains2.fasta".format(batch_index)
    FileIO.write_seqids_to_fasta(remains, cur_file, fasta_d)
    o = ar.run_minimap(cur_file, cur_file)
    print >> sys.stderr, "processing", o
    pCS, remains = sp.process_self_align_into_seed(o, remains, MiniReader, pCS)
    print >> sys.stderr, "pCS: {0}, tucked: {1}, orphans: {2}, remains: {3}".format( \
        len(pCS.S), sum(v == 'T' for v in pCS.seq_stat.itervalues()), len(orphans), len(remains))
    # combine remains+orphans to new orphans
    orphans = orphans.union(remains)
    FileIO.write_preClusterSet_to_fasta(pCS, "seed{0}.S.fasta".format(batch_index+1), fasta_d)
    FileIO.write_seqids_to_fasta(orphans, "seed{0}.orphans.fasta".format(batch_index+1), fasta_d)

    return pCS, orphans

def main():

    d = LazyFastaReader('isoseq_flnc.fasta')
    # step1. run minimap of seed0 against itself and process
    o = ar.run_minimap('seed0.fasta', 'seed0.fasta')
    seqids = set([r.id for r in SeqIO.parse(open('seed0.fasta'),'fasta')])
    pCS, orphans = sp.process_self_align_into_seed(o, seqids, MiniReader)
    # keep stats
    size_S, size_tucked, size_orphans = len(pCS.S), sum(v=='T' for v in pCS.seq_stat.itervalues()), len(orphans)
    print "seed 0 initial:"
    print size_S, size_tucked, size_orphans
    # write out seed1.S.fasta and seed1.orphans.fasta
    FileIO.write_preClusterSet_to_fasta(pCS, 'seed1.S.fasta', d)
    FileIO.write_seqids_to_fasta(orphans, 'seed1.orphans.fasta', d)
    # step 2a. minimap batch1 against seed1.S and process

    num_batchs = 5 # Liz: this is currently hardcoded for current data. CHANGE DYNAMIC LATER!!!
    for i in xrange(1, num_batchs):
        pCS, orphans = add_batch(i, pCS, orphans, d)


    # detect PCR chimeras from orphans
    chimeras = detect_PCR_chimeras(orphans, d)
    orphans = orphans.difference(chimeras)

    FileIO.write_seqids_to_fasta(orphans, "preCluster_out.orphans.fasta", d)
    FileIO.write_seqids_to_fasta(chimeras, "preCluster_out.chimeras.fasta", d)

    # dump pCS, orphans, chimeras to a pickle
    # can't dupm yet --- since pCS is an object
    #with open('preCluster.output.pickle', 'w') as f:
    #    dump({'pCS': pCS, 'chimeras': chimeras, 'orphans': orphans}, f)
    # write CSV file
    with open('preCluster.output.csv', 'w') as f:
        f.write("seqid,stat\n")
        for x, stat in pCS.seq_stat.iteritems():
            if stat == 'T': f.write("{0},tucked\n".format(x))
            elif stat == 'M': f.write("{0},{1}\n".format(x, pCS.seq_map[x]))
        for x in orphans: f.write("{0},orphan\n".format(x))
        for x in chimeras: f.write("{0},chimera\n".format(x))

    # write out a directory per preCluster cid in preCluster_out/<cid>
    singlef = open("preCluster_out.singles.fasta", 'w')
    for cid in pCS.S:
        if pCS.S[cid].size == 1:
            r = d[pCS.S[cid].members[0]]
            singlef.write(">{0}\n{1}\n".format(r.id, r.seq))
        else:
            dirname = os.path.join("preCluster_out", str(cid))
            os.makedirs(dirname)
            file = os.path.join(dirname, 'isoseq_flnc.fasta')
            FileIO.write_seqids_to_fasta(pCS.S[cid].members, file, d)
    singlef.close()

