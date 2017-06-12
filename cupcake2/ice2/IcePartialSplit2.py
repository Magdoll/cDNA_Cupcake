__author__ = 'etseng@pacb.com'

#!/usr/bin/env python
###############################################################################
# Copyright (c) 2011-2013, Pacific Biosciences of California, Inc.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
# * Neither the name of Pacific Biosciences nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
#
# NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE GRANTED BY
# THIS LICENSE.  THIS SOFTWARE IS PROVIDED BY PACIFIC BIOSCIENCES AND ITS
# CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT
# NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR
# ITS CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
# OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
# OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
# ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
###############################################################################

"""
Overview:
    pbtranscript cluster contains two main components:
    * (1) ICE (iterative clustering and error correction) to predict
      unpolished consensus isoforms.
    * (2) Polish, to use nfl reads and quiver to polish those predicted
      unpolished isoforms. Polish contains three steps:
      + (2.1) IceAllPartials (ice_partial.py all)
              Align and assign nfl reads to unpolished isoforms, and
              save results to a pickle file.
      + (2.2) IceQuiver (ice_quiver.py all)
              Call quiver to polish each isoform based on alignments
              created by mapping its associated fl and nfl reads to
              this isoform.
      + (2.3) IceQuiverPostprocess (ice_quiver.py postprocess)
              Collect and post process quiver results, and classify
              HQ/LQ isoforms.

    In order to handle subtasks by SMRTPipe instead of pbtranscript
    itself, we will refactor the polish phase including
    (2.1) (2.2) and (2.3).

    (2.1) IceAllPartials (ice_partial.py all) will be refactored to
      + (2.1.1) ice_partial.py split
                Split nfl reads into N chunks (N<=100).
      + (2.1.2) ice_partial.py i
                For each chunk i, align and assign its reads to unpolished
                isoforms and create a pickle file.
      + (2.1.3) ice_partial.py merge
                Merge pickles for all splitted chunks together to a
                big pickle.

    *** Here we are focusing on (2.1.1) ice_partial.py split ***

Description:
    (2.1.1) ice_partial.py split

    Assumption:
      * Phase (1) ICE is done and unpolished isoforms created.
      * We will split nfl reads into N chunks, N <= 100.

    Process:
        Given nfl_fa, root_dir, and N, split nfl reads into N chunks.

        The nfl reads fasta file is usually isoseq_nfl.fasta

        The i-th chunk will be placed at:
            root_dir/output/map_noFL/input.split_{0:03d}.fasta

        This job should be equivalent to command:
        fasta_splitter.py nfl_fa ceil(num_nfl_reads/N) \
                          root_dir/output/map_noFL/ \
                          input.split
    Input:
        Positional:
            nfl_fa, a fasta of non-full-length reads.
            root_dir, an output directory for running pbtranscript cluster.
            N, the number of chunks to split nfl_fa.

    Output:
        Split nfl_fa into N chunks and save.

    Hierarchy:
        pbtranscript = iceiterative

        pbtranscript --quiver = iceiterative + \
                                ice_polish.py

        ice_polish.py =  ice_make_fasta_fofn.py + \
                         ice_partial.py all + \
                         ice_quiver.py all

        ice_partial.py all = ice_partial.py split + \
                             ice_partial.py i + \
                             ice_partial.py merge

        (ice_partial.py one --> only apply ice_partial on a given input fasta)

        ice_quiver.py all = ice_quiver.py i + \
                            ice_quiver.py merge + \
                            ice_quiver.py postprocess

    Example:
        ice_partial_split.py nfl_fa root_dir N

Alternative way to call this script:
    python -m pbtranscript.ice_partial split
"""

import logging
from math import ceil
from pbtranscript.__init__ import get_version


from pbtranscript.ice.IceUtils import num_reads_in_fasta
from pbtranscript.Utils import touch, real_ppath, mkdir, nfs_exists

from cupcake2.ice2.IceFiles2 import IceFiles2
from cupcake2.io.SeqSplitter import splitFaFq
from cupcake2.tofu2.ToFuOptions2 import add_nfl_fa_argument, add_nfl_fq_argument, \
    add_cluster_root_dir_as_positional_argument
from cupcake2.ice2.__init__ import ICE_PARTIAL_PY


def add_ice_partial_split_arguments(parser):
    """Add arguments for ice_partial split."""
    parser = add_cluster_root_dir_as_positional_argument(parser)

    parser = add_nfl_fa_argument(parser, positional=True)
    parser = add_nfl_fq_argument(parser, positional=True)

    helpstr = "Split non-full-length reads into N chunks, N < 100."
    parser.add_argument("N", help=helpstr, type=int)


class IcePartialSplit2(object):

    """ice_partial split runner."""

    desc = "Split input non-full-length reads into N (N < 100) chunks, " + \
           "so that the chunked reads can be assigned later " + \
           "using '%s i' in paraellel." % ICE_PARTIAL_PY

    prog = "%s split " % ICE_PARTIAL_PY # used by cmd_str and ice_partial.py

    def __init__(self, root_dir, nfl_fa, nfl_fq, N):
        self.root_dir = root_dir
        self.nfl_fa = nfl_fa
        self.nfl_fq = nfl_fq
        self.N = int(N)
        if N >= 100:
            raise ValueError("N = {N} > 100 is too large.".format(N=self.N))

    def getVersion(self):
        """Return version string."""
        return get_version()

    def cmd_str(self):
        """Return a cmd string"""
        return self._cmd_str(root_dir=self.root_dir, nfl_fa=self.nfl_fa,
                             nfl_fq=self.nfl_fq,
                             N=self.N)

    def _cmd_str(self, root_dir, nfl_fa, nfl_fq, N):
        """Return a cmd string of ice_partial_split."""
        cmd = self.prog + \
            "{d} ".format(d=root_dir) + \
            "{nfl_fa} ".format(nfl_fa=nfl_fa) + \
            "{nfl_fq} ".format(nfl_fq=nfl_fq) + \
            "{N} ".format(N=N)
        return cmd

    def validate_inputs(self):
        """
        Check whether root_dir exists, split nfl reads in nfl_fa into N chunks.
        """
        return self._validate_inputs(root_dir=self.root_dir, nfl_fa=self.nfl_fa,
                                     nfl_fq=self.nfl_fq,
                                     N=self.N)

    def _validate_inputs(self, root_dir, nfl_fa, nfl_fq, N):
        """
        Check inputs, return
        (num_reads,
         number_reads_per_chunk,
         nfl_dir,
         [i-th_chunk_nfl_fa for i in [0...N-1]])
        """
        icef = IceFiles2(prog_name="ice_partial_split",
                        root_dir=root_dir, no_log_f=False)

        nfl_dir = icef.nfl_dir

        # root_dir/output/map_noFL/input.split_{0:03d}.fasta
        splitted_nfl_fas = [icef.nfl_fa_i(i) for i in range(0, N)]
        splitted_nfl_fqs = [icef.nfl_fq_i(i) for i in range(0, N)]

        mkdir(icef.nfl_dir)

        # Check if inputs exist.
        errMsg = ""
        if N <= 0 or N > 100:
            errMsg = "Input file can not be splitted into %d chunks!" % N

        if not nfs_exists(nfl_fa):
            errMsg = ("The input non-full-length reads fasta file " +
                      "{f} does not exists. ".format(f=nfl_fa))
        if not nfs_exists(nfl_fq):
            errMsg = ("The input non-full-length reads fastq file " +
                      "{f} does not exists. ".format(f=nfl_fq))
        if len(errMsg) != 0:
            raise ValueError(errMsg)

        num_reads = num_reads_in_fasta(nfl_fa)
        reads_per_split = int(max(1, ceil(num_reads * 1.0 / N)))

        return (num_reads, reads_per_split, nfl_dir, splitted_nfl_fas, splitted_nfl_fqs)

    def run(self):
        """Run"""
        logging.debug("root_dir: {d}.".format(d=self.root_dir))
        logging.debug("nfl_fa: {f}.".format(f=self.nfl_fa))
        logging.debug("nfl_fq: {f}.".format(f=self.nfl_fq))
        logging.debug("Total number of chunks: N={N}.".format(N=self.N))

        # Validate input files,
        (num_reads, reads_per_split, nfl_dir, splitted_fas_todo, splitted_fqs_todo) = \
            self.validate_inputs()

        logging.info("Total number of reads is {n}.".format(n=num_reads))
        logging.info("Splitting nfl_fa into chunks each " +
                     "containing {n} reads.".format(n=reads_per_split))

        splitted_fas_done = splitFaFq(
            input_fa_or_fq=self.nfl_fa,
            reads_per_split=reads_per_split,
            out_dir=nfl_dir,
            out_format=IceFiles2.nfl_fa_format,
            is_fq=False,
            reads_in_first_split=None)

        splitted_fqs_done = splitFaFq(
            input_fa_or_fq=self.nfl_fq,
            reads_per_split=reads_per_split,
            out_dir=nfl_dir,
            out_format=IceFiles2.nfl_fq_format,
            is_fq=True,
            reads_in_first_split=None)

        logging.info("Splitted fastas are: " + "\n".join(splitted_fas_done))
        logging.info("Splitted fastqs are: " + "\n".join(splitted_fqs_done))

        # Liz: why is below necessary????
        for fa in splitted_fas_todo:
            if fa not in splitted_fas_done:
                logging.info("touching {f}".format(f=fa))
                touch(fa)
        for fq in splitted_fqs_todo:
            if fq not in splitted_fqs_done:
                logging.info("touching {f}".format(f=fq))
                touch(fq)
