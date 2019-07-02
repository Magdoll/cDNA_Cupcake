__author__ = 'etseng@pacb.com'


import os, sys, pdb, subprocess
from collections import defaultdict
from functools import partial

import networkx as nx
from Bio import SeqIO

from cupcake2.ice2.preCluster import preClusterSet2

def sanity_checking(pCS, orphans):
    # do some sanity checking
    # 1. all members should be uniquely in one S
    for cid in pCS.S:
        for x in pCS.S[cid].members:
            assert pCS.seq_map[x]==cid
            assert x not in orphans
    # 2. anything in tucked must not be in seq_map
    for x in pCS.tucked:
        assert x not in pCS.seq_map
        assert x not in orphans

def sanity_checking2(pCS, orphans):
    for cid in pCS.S:
        for x in pCS.S[cid].members:
            assert pCS.seq_map[x]==cid
            assert x not in orphans

    for x in pCS.seq_map:
        if pCS.seq_stat[x] == 'T':
            assert x not in orphans
            assert x not in pCS.seq_map
        elif pCS.seq_stat[x] == 'M':
            assert x not in orphans


def process_self_align_into_seed(align_filename, seqids, reader_class, pCS=None, dun_use_partial=False):
    """
    Ignore hits that are - strand or qID >= sID (self hit or already reported)

    Returns:
    pCS -- preClusterSet
    orphans --- seqs that are neither in pCS nor tucked i.e. no align hits
    """
    if pCS is None:
        pCS = preClusterSet2()

    orphans = set(seqids)
    reader = reader_class(align_filename)
    for r in reader:
        if r.qID >= r.sID or r.strand == '-': continue
        s = r.characterize(30, 0.01, 30, 0.01, 30, 0.05, min_identity=0.99)
        if dun_use_partial and s == 'partial': continue

        if s == 'match':
            pCS.add_seqid_match(r.qID, r.sID)
            if r.qID in bug_ids or r.sID in bug_ids:
                print "match:", r.qID, r.sID
                print "after match:", pCS.seq_stat[r.qID], pCS.seq_stat[r.sID]
        elif s == 'partial':
            pCS.add_seqid_partial(r.qID, r.sID)
            if r.qID in bug_ids or r.sID in bug_ids:
                print "partial:", r.qID, r.sID
                print "after partial:", pCS.seq_stat[r.qID], pCS.seq_stat[r.sID]
        elif s == 'q_contained':
            if r.qID in bug_ids or r.sID in bug_ids:
                print "tucking {0} into {1}".format(r.qID, r.sID)
                print "before:", pCS.seq_stat[r.qID], pCS.seq_stat[r.sID]
            pCS.add_seqid_contained(r.qID, r.sID)
            if r.qID in bug_ids or r.sID in bug_ids:
                print "after:", pCS.seq_stat[r.qID], pCS.seq_stat[r.sID]
        elif s == 's_contained':
            if r.qID in bug_ids or r.sID in bug_ids:
                print "tucking {0} into {1}".format(r.sID, r.qID)
                print "before:", pCS.seq_stat[r.sID], pCS.seq_stat[r.qID]
            pCS.add_seqid_contained(r.sID, r.qID)
            if r.qID in bug_ids or r.sID in bug_ids:
                print "after:", pCS.seq_stat[r.sID], pCS.seq_stat[r.qID]
        try:
            orphans.remove(r.qID)
        except:
            pass
        try:
            orphans.remove(r.sID)
        except:
            pass
        #sanity_checking(pCS, orphans)


    #sanity_checking(pCS, orphans)

    return pCS, orphans

def process_align_to_pCS(align_filename, seqids, pCS, reader_class, dun_use_partial=False):
    """
    Batch against {S}, so it is NOT self align! Only have to ignore - strand hits.

    (a) if match against exactly one cluster c: add to c
    (b) if match against multiple clusters: combine (add_seqid_match handles this)
    (c) if match against none: save as "batch orphan" to be processed in step 2b

    Return: pCS, remains
    """
    orphans = set(seqids)
    # ex: process batch1 against seed1.S.fasta
    reader = reader_class(align_filename)#'batch1.fasta.S.f00001.minimap')
    for r in reader:
        if r.strand == '-': continue
        s = r.characterize(30, 0.01, 30, 0.01, 30, 0.05, min_identity=0.99)
        if dun_use_partial and s == 'partial': continue
        # Liz note: currently, just add all to match because minimap sensitivity not enough to do "tuck" properly
        if s == 'match':
            pCS.add_seqid_match(r.qID, r.sID)
            if r.qID in bug_ids or r.sID in bug_ids:
                print "match:", r.qID, r.sID
                print "after match:", pCS.seq_stat[r.qID], pCS.seq_stat[r.sID]
        elif s == 'partial':
            if r.qID in bug_ids or r.sID in bug_ids:
                print "partial:", r.qID, r.sID
                print "after partial:", pCS.seq_stat[r.qID], pCS.seq_stat[r.sID]
            pCS.add_seqid_partial(r.qID, r.sID)
        elif s == 'q_contained':
            # sID must be in cluster, so just call pCS to handle the tucking
            if r.qID in bug_ids or r.sID in bug_ids:
                print "tucking {0} into {1}".format(r.qID, r.sID)
                print "before:", pCS.seq_stat[r.qID], pCS.seq_stat[r.sID]
            pCS.add_seqid_contained(r.qID, r.sID)
            if r.qID in bug_ids or r.sID in bug_ids:
                print "after:", pCS.seq_stat[r.qID], pCS.seq_stat[r.sID]
        elif s == 's_contained':
            # sID must be in cluster
            if r.qID in bug_ids or r.sID in bug_ids:
                print "tucking {0} into {1}".format(r.sID, r.qID)
                print "before:", pCS.seq_stat[r.sID], pCS.seq_stat[r.qID]
            pCS.add_seqid_contained(r.sID, r.qID)
            if r.qID in bug_ids or r.sID in bug_ids:
                print "after:", pCS.seq_stat[r.sID], pCS.seq_stat[r.qID]
        try:
            orphans.remove(r.qID)
        except:
            pass

    return pCS, orphans


bug_ids = [line.strip() for line in open('ids_for_NitrateReductase_gene.txt')]

def process_align_to_orphan(align_filename, remaining, orphans, pCS, reader_class, dun_use_partial=False):
    """
    remaining of "batch" against orphans

    if partial: ignore
    if match: add orphan+matched to {S} (remember orphan is not currently in S); also remove orphan from "orphan"; remove qID from tmp
    if q_contained: add orphan to S, tuck qID; also remove orphan from the variable "orphan"; remove qID from tmp
    if s_contained: add qID to S, tuck orphan; also remove orphan from the variable "orphan"; remove qID from tmp

    Returns: pCS, tucked, orphans, remaining
    """
    reader = reader_class(align_filename)
    for r in reader:
        if r.strand == '-': continue
        s = r.characterize(30, 0.01, 30, 0.01, 30, 0.05, min_identity=0.99)
        if dun_use_partial and s == 'partial': continue
        if s == 'match':
            pCS.add_seqid_match(r.qID, r.sID)
            if r.qID in bug_ids or r.sID in bug_ids:
                print "after match:", pCS.seq_stat[r.qID], pCS.seq_stat[r.sID]
        elif s == 'partial':
            pCS.add_seqid_partial(r.qID, r.sID)
            if r.qID in bug_ids or r.sID in bug_ids:
                print "after partial:", pCS.seq_stat[r.qID], pCS.seq_stat[r.sID]
        elif s == 'q_contained':
            if r.qID in bug_ids or r.sID in bug_ids:
                print "tucking {0} into {1}".format(r.qID, r.sID)
                print "before:", pCS.seq_stat[r.qID], pCS.seq_stat[r.sID]
            pCS.add_seqid_contained(r.qID, r.sID)
            if r.qID in bug_ids or r.sID in bug_ids:
                print "after:", pCS.seq_stat[r.qID], pCS.seq_stat[r.sID]
        elif s == 's_contained':
            if r.qID in bug_ids or r.sID in bug_ids:
                print "tucking {0} into {1}".format(r.sID, r.qID)
                print "before:", pCS.seq_stat[r.sID], pCS.seq_stat[r.qID]
            pCS.add_seqid_contained(r.sID, r.qID)
            if r.qID in bug_ids or r.sID in bug_ids:
                print "after:", pCS.seq_stat[r.sID], pCS.seq_stat[r.qID]
        try:
            orphans.remove(r.sID)
        except:
            pass

        try:
            remaining.remove(r.qID)
        except:
            pass

    return pCS, orphans, remaining