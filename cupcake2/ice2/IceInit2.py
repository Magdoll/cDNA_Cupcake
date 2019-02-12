#!/usr/bin/env python
"""
Find initial mutual exclusive cliques by aligning
input reads against itself.

This is for ICE2 --- as of 2019, it uses ONLY minimap2.


Several differences with the old IceInit:

1. Only uses minimap2 (no more BLASR or DALIGNER)!
2. All QV usages are removed
3. Dependencies to old pbtranscript/pbtranscript2 removed.

"""
__author__ = 'etseng@pacb.com'


import time
import logging
import networkx as nx
from Bio import SeqIO

import cupcake.ice.pClique as pClique
from cupcake2.ice2.AlignerRunners import run_minimap
from cupcake2.ice2.IceUtils2 import minimap2_against_ref2

class IceInit2(object):
    """Iterative clustering and error correction."""
    def __init__(self, readsFa,
                 ice_opts, sge_opts, run_uc=True):

        self.readsFa = readsFa
        self.ice_opts = ice_opts
        self.sge_opts = sge_opts


        self.ice_opts.detect_cDNA_size(readsFa)

        self.uc = None
        if run_uc:
            self.uc = self.init_cluster_by_clique()


    @classmethod
    def _findCliques(self, alignGraph, readsFa):
        """
        Find all mutually exclusive cliques within the graph, with decreased
        size.

        alignGraph - a graph, each node represent a read and each edge
        represents an alignment between two end points.

        Return a dictionary of clique indices and nodes.
            key = index of a clique
            value = nodes within a clique
        Cliques are ordered by their size descendingly: index up, size down
        Reads which are not included in any cliques will be added as cliques
        of size 1.
        """
        uc = {}    # To keep cliques found
        used = []  # nodes within any cliques
        ind = 0    # index of clique to discover

        deg = list(alignGraph.degree())
        # Sort tuples of (node, degree) by degree, descendingly
        deg.sort(key=lambda x: x[1], reverse=True)
        for d in deg:
            node = d[0]  # node which has the largest degree in alignGraph
            if node not in alignGraph:
                continue
            # just get the immediate neighbors since we're looking for perfect
            # cliques
            subGraph = alignGraph.subgraph([node] + list(alignGraph.neighbors(node)))
            subNodes = list(subGraph.nodes())
            # Convert from networkx.Graph to a sparse matrix
            S, H = pClique.convert_graph_connectivity_to_sparse(
                subGraph, subNodes)
            # index of the 'node' in the sub-graph
            seed_i = subNodes.index(node)
            # Grasp a clique from subGraph, and return indices of clique nodes
            # setting gamma=0.8 means to find quasi-0.8-cliques!
            tQ = pClique.grasp(S, H, gamma=0.8, maxitr=5, given_starting_node=seed_i)
            if len(tQ) > 0:
                c = [subNodes[i] for i in tQ]  # nodes in the clique
                uc[ind] = c  # Add the clique to uc
                ind += 1
                used += c    # Add clique nodes to used
                # Remove clique nodes from alignGraph and continue
                alignGraph.remove_nodes_from(c)

        # write each orphan as a singleton cluster
        for r in SeqIO.parse(open(readsFa), 'fasta'):
            if r.id not in used:
                uc[ind] = [r.id]
                ind += 1
        return uc

    def makeGraphFromMinimap2(self, align_filename, len_dict):
        alignGraph = nx.Graph()

        edge_count = 0
        start_t = time.time()
        for r in minimap2_against_ref2(
                sam_filename=align_filename,
                query_len_dict=len_dict,
                ref_len_dict=len_dict,
                is_FL=True,
                sID_starts_with_c=False,
                ece_penalty=self.ice_opts.ece_penalty,
                ece_min_len=self.ice_opts.ece_min_len,
                same_strand_only=True,
                max_missed_start=self.ice_opts.max_missed_start,
                max_missed_end=self.ice_opts.max_missed_end,
                full_missed_start=self.ice_opts.full_missed_start,
                full_missed_end=self.ice_opts.full_missed_end):
            if r.qID == r.cID:
                continue  # self hit, ignore
            if r.ece_arr is not None:
                alignGraph.add_edge(r.qID, r.cID)
                edge_count += 1

        logging.debug("total {0} edges added from {1}; took {2} sec".format(\
            edge_count, align_filename, time.time() - start_t))
        return alignGraph



    def init_cluster_by_clique(self):
        """
        Only called once and in the very beginning, when (probably a subset)
        of sequences are given to generate the initial cluster.

        Returns dict of cluster_index --> list of seqids
        which is the 'uc' dict that can be used by IceIterative
        """

        # run minimap2 of self vs self
        align_filename = run_minimap(self.readsFa, self.readsFa, \
                                     cpus=self.sge_opts.blasr_nproc, \
                                     sam_output=True)
        len_dict = dict((r.id, len(r.seq)) for r in SeqIO.parse(open(self.readsFa), 'fasta'))
        alignGraph = self.makeGraphFromMinimap2(align_filename, len_dict)

        uc = IceInit2._findCliques(alignGraph=alignGraph, readsFa=self.readsFa)
        return uc

