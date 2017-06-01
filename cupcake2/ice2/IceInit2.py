#!/usr/bin/env python
"""
Find initial mutual exclusive cliques by aligning
input reads against itself.

This is for ICE2 --- it uses either DALIGNER (default) or BLASR.


Several differences with the old IceInit:

1. It can be run in a standalone manner! Does not have to tie with the whole ICE2 framework.
2. It does NOT use pbcore. All pbcore is replaced with BioPython or other functions.
3. The BLASR version is expected to be SA4.0+, so options are -- not -.
4. It does not use QVs!

"""
__author__ = 'etseng@pacb.com'

import os
import os.path as op
import time
import logging
import subprocess
import networkx as nx
import numpy as np
from Bio import SeqIO

import pbtranscript.ice.pClique as pClique
from pbtranscript.ice.IceUtils import blasr_against_ref, daligner_against_ref
from pbtranscript.ice.IceUtils import set_probqv_from_model
from pbtranscript.Utils import real_upath, execute
from pbtranscript.ice_daligner import DalignerRunner


class IceInit2(object):
    """Iterative clustering and error correction."""
    def __init__(self, readsFa, qver_get_func, qvmean_get_func,
                 ice_opts, sge_opts, run_uc=True):

        self.readsFa = readsFa
        self.ice_opts = ice_opts
        self.sge_opts = sge_opts

        self.qver_get_func = qver_get_func
        self.qvmean_get_func = qvmean_get_func

        self.ice_opts.detect_cDNA_size(readsFa)

        self.uc = None
        if run_uc:
            self.uc = self.init_cluster_by_clique()

    # version using BLASR; fallback if daligner fails
    def _align_withBLASR(self, queryFa, targetFa, outFN):
        """Align input reads against itself using BLASR."""
        if op.exists(outFN):
            logging.info("{0} already exists. No need to run BLASR.".format(outFN))
        else:
            cmd = "blasr {q} ".format(q=real_upath(queryFa)) + \
                  "{t} ".format(t=real_upath(targetFa)) + \
                  "-m 5 --maxLCPLength 15 " + \
                  "--nproc {cpu} ".format(cpu=self.sge_opts.blasr_nproc) + \
                  "--minAlnLength {aln} ".format(aln=self.ice_opts.min_match_len) + \
                  "--maxScore {score} ".format(score=self.ice_opts.maxScore) + \
                  "--bestn {n} --nCandidates {n} ".format(n=self.ice_opts.bestn) + \
                  "--out {o} ".format(o=real_upath(outFN)) + \
                  "1>/dev/null 2>/dev/null"
            logging.info("Calling {cmd}".format(cmd=cmd))
            execute(cmd)

    # align with DALIGNER
    def _align_withDALIGNER(self, queryFa, output_dir):
        """
        Align input reads against itself using DALIGNER.
        """
        # run this locally
        # Liz: is_FL is currently turned OFF! because LA4Ice has ICE_FL(-E) set with 200/50bp missed, too strict
        runner = DalignerRunner(query_filename=queryFa, target_filename=queryFa,
                                query_converted=False, target_converted=False,
                                is_FL=False, same_strand_only=True,
                                use_sge=False, sge_opts=None,
                                cpus=4)

        runner.run(min_match_len=self.ice_opts.min_match_len,
                   output_dir=output_dir,
                   sensitive_mode=self.ice_opts.sensitive_mode)
        return runner

    # version using BLASR
    def _makeGraphFromM5(self, m5FN):
        """Construct a graph from a BLASR M5 file."""
        alignGraph = nx.Graph()

        for r in blasr_against_ref(output_filename=m5FN,
            is_FL=True,
            sID_starts_with_c=False,
            qver_get_func=self.qver_get_func,
            qvmean_get_func=self.qvmean_get_func,
            ece_penalty=self.ice_opts.ece_penalty,
            ece_min_len=self.ice_opts.ece_min_len,
            max_missed_start=self.ice_opts.max_missed_start,
            max_missed_end=self.ice_opts.max_missed_end):
            if r.qID == r.cID:
                continue # self hit, ignore
            if r.ece_arr is not None:
                logging.debug("adding edge {0},{1}".format(r.qID, r.cID))
                alignGraph.add_edge(r.qID, r.cID)
        return alignGraph


    def _makeGraphFromLA4Ice(self, runner):
        """Construct a graph from a LA4Ice output file."""
        alignGraph = nx.Graph()

        for la4ice_filename in runner.la4ice_filenames:
            count = 0
            start_t = time.time()
            for r in daligner_against_ref(
                    query_dazz_handler=runner.query_dazz_handler,
                    target_dazz_handler=runner.target_dazz_handler,
                    la4ice_filename=la4ice_filename,
                    is_FL=True, sID_starts_with_c=False,
                    qver_get_func=self.qver_get_func, qvmean_get_func=self.qvmean_get_func,
                    qv_prob_threshold=.03, ece_min_len=self.ice_opts.ece_min_len,
                    ece_penalty=self.ice_opts.ece_penalty,
                    same_strand_only=True, no_qv_or_aln_checking=False,
                    max_missed_start=self.ice_opts.max_missed_start,
                    max_missed_end=self.ice_opts.max_missed_end):
                if r.qID == r.cID:
                    continue # self hit, ignore
                if r.ece_arr is not None:
                    alignGraph.add_edge(r.qID, r.cID)
                    count += 1
            logging.debug("total {0} edges added from {1}; took {2} sec"
                          .format(count, la4ice_filename, time.time()-start_t))
        return alignGraph


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

        deg = alignGraph.degree().items()
        # Sort tuples of (node, degree) by degree, descendingly
        deg.sort(key=lambda x: x[1], reverse=True)
        for d in deg:
            node = d[0]  # node which has the largest degree in alignGraph
            if node not in alignGraph:
                continue
            # just get the immediate neighbors since we're looking for perfect
            # cliques
            subGraph = alignGraph.subgraph([node] + alignGraph.neighbors(node))
            subNodes = subGraph.nodes()
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


    def init_cluster_by_clique(self):
        """
        Only called once and in the very beginning, when (probably a subset)
        of sequences are given to generate the initial cluster.

        readsFa --- initial fasta filename, probably called *_split00.fasta
        qver_get_func --- function that returns QVs on reads
        qvmean_get_func --- function that returns the mean QV on reads
        bestn --- parameter in BLASR, higher helps in finding perfect
            cliques but bigger output
        nproc, maxScore --- parameter in BLASR, set maxScore appropriate
            to input transcript length
        ece_penalty, ece_min_len --- parameter in isoform hit calling

        Self-blasr input then iteratively find all mutually exclusive
            cliques (in decreasing size)
        Returns dict of cluster_index --> list of seqids
        which is the 'uc' dict that can be used by IceIterative
        """
        alignGraph = None

        if self.ice_opts.aligner_choice == 'blasr':
            outFN = self.readsFa + '.self.blasr'
            self._align_withBLASR(queryFa=self.readsFa, targetFa=self.readsFa, outFN=outFN)
            alignGraph = self._makeGraphFromM5(m5FN=outFN)
        elif self.ice_opts.aligner_choice == 'daligner':
            try:
                runner = self._align_withDALIGNER(queryFa=self.readsFa,
                                                  output_dir=op.dirname(real_upath(self.readsFa)))
                alignGraph = self._makeGraphFromLA4Ice(runner=runner)
                runner.clean_run()
            except RuntimeError:  # daligner probably crashed, fall back to blasr
                outFN = self.readsFa + '.self.blasr'
                self._align_withBLASR(queryFa=self.readsFa, targetFa=self.readsFa, outFN=outFN)
                alignGraph = self._makeGraphFromM5(m5FN=outFN)
        else:
            raise Exception, "Unrecognized aligner_choice {0}!".format(self.ice_opts.aligner_choice)

        uc = IceInit2._findCliques(alignGraph=alignGraph, readsFa=self.readsFa)
        return uc

