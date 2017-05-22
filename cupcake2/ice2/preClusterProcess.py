__author__ = 'etseng@pacb.com'


import os, sys, pdb, subprocess
from collections import defaultdict

import networkx as nx
from Bio import SeqIO

from cupcake2.ice2.preCluster import preClusterSet
from cupcake2.io.minimapIO import MiniReader


def process_self_align_into_pCS(align_filename, seqids, reader_class=MiniReader):
    """
    Ignore hits that are - strand or qID >= sID (self hit or already reported)

    Returns:
    pCS -- preClusterSet
    tucked --- dict of seqid --> what it is tucked under
    orphans --- seqs that are neither in pCS nor tucked i.e. no align hits
    """
    pCS = preClusterSet()
    tucked = defaultdict(lambda: []) # seqid --> what it is tucked under
    orphans = set(seqids)
    reader = reader_class(align_filename)
    for r in reader:
        if r.qID >= r.sID or r.strand == '-': continue
        s = r.characterize(400, 0.1, 400, 0.1)
        if s == 'partial': continue
        if s == 'match': pCS.add_seqid_match(r.qID, r.sID)
        elif s == 'q_contained': tucked[r.qID].append(r.sID)
        elif s == 's_contained': tucked[r.sID].append(r.qID)

    # use networkx to get the "tucked" map ones that have no outgoing nodes (but also mean they have no matches)
    G = nx.DiGraph()
    for k,v in tucked.iteritems():
        for s in v:
            G.add_edge(k, s)

    cand = filter(lambda x: G.out_degree(x)==0 and x not in pCS.seq_map, G.nodes_iter())
    # add candidates (from tucked, no outgoing, not in pCS already) to pCS.S
    for x in cand: pCS.add_new_cluster([x])


    # do some sanity checking
    # 1. all members should be uniquely in one S
    for cid in pCS.S:
        for x in pCS.S[cid].members:
            assert pCS.seq_map[x]==cid
    # 2. get the set of "orphans" (not in tucked, not in pCS.S)
    for x in pCS.seq_map: orphans.remove(x)
    for x in tucked:
        try: orphans.remove(x)
        except: pass # ignore, already removed

    return pCS, tucked, orphans