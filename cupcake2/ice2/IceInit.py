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
import time
import logging
import subprocess
import networkx as nx
import numpy as np
from Bio import SeqIO

import pbtools.pbtranscript.ice.pClique as pClique
from pbtools.pbtranscript.ice.IceUtils import blasr_against_ref, get_daligner_sensitivity_setting
from pbtools.pbtranscript.icedalign.IceDalignUtils import DazzIDHandler, DalignerRunner
from pbtools.pbtranscript.icedalign.IceDalignReader import dalign_against_ref

class IceInit(object):
    """Iterative clustering and error correction."""
    def __init__(self, readsFa, ice_opts, sge_opts):

        self.readsFa = readsFa
        self.ice_opts = ice_opts
        self.sge_opts = sge_opts
        self.uc = None

        self.uc = self.init_cluster_by_clique(
            readsFa=readsFa,
            ice_opts=ice_opts, sge_opts=sge_opts)


    # version using BLASR; fallback if daligner fails
    def _align_withBLASR(self, queryFa, targetFa, outFN, ice_opts, sge_opts):
        """Align input reads against itself using BLASR."""
        if os.path.exists(outFN):
            logging.info("{0} already exists. No need to run BLASR.".
                format(outFN))
        else:
            cmd = "blasr {q} ".format(q=queryFa) + \
                  "{t} ".format(t=targetFa) + \
                  "-m 5 --maxLCPLength 15 " + \
                  "--nproc {cpu} ".format(cpu=sge_opts.blasr_nproc) + \
                  "--maxScore {score} ".format(score=ice_opts.maxScore) + \
                  "--bestn {n1} --nCandidates {n2} ".format(n1=ice_opts.bestn, n2=ice_opts.bestn*2) + \
                  "-out {o}".format(o=outFN)
            logging.info("Calling {cmd}".format(cmd=cmd))

            try:
                code = subprocess.check_call(cmd, shell=True)
                if code!=0:
                    errMsg = "{cmd} exited with {code}".\
                        format(cmd=cmd, code=code)
                    logging.error(errMsg)
                    raise RuntimeError(errMsg)
            except subprocess.CalledProcessError as grepexc:
                errMsg = "{cmd} exited with {code}:{msg}".format(cmd=cmd,\
                    code=grepexc.returncode, msg=grepexc.output)
                raise RuntimeError(errMsg)

    # using DALIGNER, fallback is _align_with_BLASR
    def _align(self, queryFa, output_dir, ice_opts, sge_opts):

        daligner_sensitive_mode, _low, _high, _ignore5, _ignore3, _ece_min_len, _ece_penalty = \
            get_daligner_sensitivity_setting(queryFa, is_fasta=True)

        input_obj = DazzIDHandler(queryFa, False)
        DalignerRunner.make_db(input_obj.dazz_filename)

        # run this locally
        runner = DalignerRunner(queryFa, queryFa, is_FL=True, same_strand_only=True, \
                            query_converted=True, db_converted=True, query_made=True, \
                            db_made=True, use_sge=False, cpus=4, sge_opts=None)
        las_filenames, las_out_filenames = runner.runHPC(min_match_len=_low, output_dir=output_dir, sensitive_mode=daligner_sensitive_mode)
        return input_obj, las_out_filenames, _ignore5, _ignore3, _ece_min_len, _ece_penalty


    # version using BLASR
    def _makeGraphFromM5(self, m5FN, qver_get_func, qvmean_get_func, ice_opts, max_missed_start, max_missed_end):
        """Construct a graph from a BLASR M5 file."""
        alignGraph = nx.Graph()

        for r in blasr_against_ref(output_filename=m5FN,
            is_FL=True,
            sID_starts_with_c=False,
            qver_get_func=qver_get_func,
            qvmean_get_func=qvmean_get_func,
            ece_penalty=ice_opts.ece_penalty,
            ece_min_len=ice_opts.ece_min_len,
            max_missed_start=max_missed_start,
            max_missed_end=max_missed_end):
            if r.qID == r.cID:
                continue # self hit, ignore
            if r.ece_arr is not None:
                logging.debug("adding edge {0},{1}".format(r.qID, r.cID))
                alignGraph.add_edge(r.qID, r.cID)
        return alignGraph

    def _makeGraphFromLasOut(self, las_out_filenames, dazz_obj, qver_get_func, ice_opts, max_missed_start, max_missed_end, qvmean_get_func=None):
        alignGraph = nx.Graph()

        for las_out_filename in las_out_filenames:
            count = 0
            start_t = time.time()
            for r in dalign_against_ref(dazz_obj, dazz_obj, las_out_filename, is_FL=True, sID_starts_with_c=False,
                      qver_get_func=qver_get_func, qvmean_get_func=qvmean_get_func,
                      qv_prob_threshold=.03,
                      ece_penalty=ice_opts.ece_penalty, ece_min_len=ice_opts.ece_min_len,
                      same_strand_only=True, no_qv_or_aln_checking=False,
                      max_missed_start=max_missed_start, max_missed_end=max_missed_end):
                if r.qID == r.cID: continue # self hit, ignore
                if r.ece_arr is not None:
                    alignGraph.add_edge(r.qID, r.cID)
                    count += 1
            logging.debug("total {0} edges added from {1}; took {2} sec".format(count, las_out_filename, time.time()-start_t))
        return alignGraph


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
        deg.sort(key=lambda x:x[1], reverse=True)
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

        with SeqIO.parse(readsFa, 'fasta') as reader:
            for r in reader:
                rid = r.name.split()[0]
                if rid not in used:
                    uc[ind] = [rid]
                    ind += 1
        return uc

    def init_cluster_by_clique(self, readsFa, ice_opts, sge_opts):
        """
        Only called once and in the very beginning, when (probably a subset)
        of sequences are given to generate the initial cluster.

        readsFa --- initial fasta filename, probably called *_split00.fa
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
        output_dir = os.path.dirname(readsFa)

        try:
            dazz_obj, las_out_filenames, _ignore5, _ignore3, _ece_min_len = self._align(queryFa=readsFa, output_dir=output_dir,
                ice_opts=ice_opts, sge_opts=sge_opts)
            ice_opts.ece_min_len = _ece_min_len # overwrite
            alignGraph = self._makeGraphFromLasOut(las_out_filenames, dazz_obj, qver_get_func, ice_opts, \
                                                   max_missed_start=_ignore5, max_missed_end=_ignore3, \
                                                   qvmean_get_func=qvmean_get_func)
        except: # daligner probably crashed, fall back to blasr
            outFN = readsFa + '.self.blasr'
            daligner_sensitive_mode, _low, _high, _ignore5, _ignore3, _ece_min_len = get_daligner_sensitivity_setting(readsFa, is_fasta=True)
            ice_opts.ece_min_len = _ece_min_len # overwrite
            self._align_withBLASR(queryFa=readsFa, targetFa=readsFa, outFN=outFN, ice_opts=ice_opts, sge_opts=sge_opts)
            alignGraph = self._makeGraphFromM5(outFN, qver_get_func, qvmean_get_func, ice_opts, \
                                               max_missed_start=_ignore5, max_missed_end=_ignore3)

        uc = self._findCliques(alignGraph=alignGraph, readsFa=readsFa)
        return uc



