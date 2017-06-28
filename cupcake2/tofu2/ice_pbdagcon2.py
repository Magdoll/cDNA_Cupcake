#!/usr/bin/env python

###############################################################################
# Copyright (c) 2011-2014, Pacific Biosciences of California, Inc.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
# * Neither the name of Pacific Biosciences nor the names of its contributors
#   may be used to endorse or promote products derived from this software
#   without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY PACIFIC BIOSCIENCES AND ITS CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
# ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
###############################################################################

"""
This is a pbdagcon wrapper for pbtranscript.

Precursor to using the smrtpipe pbdagcon
(1) Find the best seed as reference
(2) Align rest to seed
(3) Call pbdagcon
"""

from argparse import ArgumentParser
from collections import defaultdict
import os
import sys

import numpy as np

from pbcore.io import FastaReader
from pbcore.util.Process import backticks
from pbtranscript.io import FastaRandomReader
from pbtranscript.__init__ import get_version

__author__ = 'etseng@pacificbiosciences.com'

class AlignGraphUtilError(Exception):
    """Align Group Util Error Class"""
    pass


def choose_template_by_blasr(fasta_filename, out_filename, nproc=8,
                             maxScore=-1000, min_number_reads=1):
    """
    Choose the best template for gcon reference
    Pick the one that has the highest average hit identity to others

    Returns: FastaRecord of selected ref
    """
    fd = FastaRandomReader(fasta_filename)

    cmd = "blasr --nproc {nproc} ".format(nproc=nproc) + \
          "--maxScore {score} ".format(score=maxScore) + \
          "--maxLCPLength 15 --bestn 10 --nCandidates 50 " + \
          "-m 1 {fa} {fa} ".format(fa=fasta_filename) + \
          "--out {out} ".format(out=out_filename) + \
          "1>/dev/null 2>/dev/null"

    out, code, msg = backticks(cmd)
    if code != 0:
        return None

    # blasr -m 1 output format:
    # (0) qName   (1) tName (2) qStrand
    # (3) tStrand (4) score (5) percentSimilarity
    # (6) tStart  (7) tEnd  (8) tLength
    # (9) qStart  (10) qEnd (11) qLength
    # (12) nCells
    scores = defaultdict(lambda: [])
    with open(out_filename) as f:
        for line in f:
            raw = line.strip().split()
            # qID gets an extra /0_length
            qID, tID = raw[0][:raw[0].rfind('/')], raw[1]
            if qID == tID:
                continue  # self-hit, ignore
            if raw[2] != raw[3]:
                continue  # has to be on same strand
            scores[qID].append(abs(float(raw[4])))  # Liz: use score as the scorer!

    # find the one with the highest average alignment similarity
    score_array = []
    for k, v in scores.iteritems():
        score_array.append((np.ceil(np.mean(v)), k))
    if len(score_array) < min_number_reads:
        errMsg = "Not enough number of reads in " + \
                 "choose_template_by_blasr {0} < {1}".format(
                     len(score_array), min_number_reads)
        raise AlignGraphUtilError(errMsg)

    score_array.sort(reverse=True)

    # Find the longest sequence that is within the std deviation of
    # the best score
    best_mean_std = np.std([x[0] for x in score_array])
    best_mean, best_id = score_array[0]
    best_len = len(fd[best_id].sequence)
    for _mean, _id in score_array[1:]:
        if abs(_mean-best_mean) > best_mean_std:
            break
        _len = len(fd[_id].sequence)
        if _len > best_len:
            best_id = _id
            best_len = _len

    return fd[best_id]


def make_aln_input_to_ref(fasta_filename, ref_filename,
                          out_filename, nproc=8):
    """
    Make blasr -m 5 output of input aligned to ref

    Return: alignment filename
    """
    # pbdagcon only takes -m 5 output
    tmp_out = "{out}.tmp".format(out=out_filename)

    cmd = "blasr {infa} ".format(infa=fasta_filename) + \
          "{ref} --bestn 1 ".format(ref=ref_filename) + \
          "--nproc {nproc} ".format(nproc=nproc) + \
          "-m 5 --out {out} ".format(out=tmp_out) + \
          "1>/dev/null 2>/dev/null"

    out, code, msg = backticks(cmd)
    if code != 0:
        errMsg = "Unable to align {infa} to {ref} in make_aln_input_to_ref".\
            format(infa=fasta_filename, ref=ref_filename)
        raise AlignGraphUtilError(errMsg)

    # trim away anything that is a self-hit or opp strand
    with open(out_filename, 'w') as f, \
            open(tmp_out, 'r') as h:
        for line in h:
            raw = line.strip().split()
            # blasr -m 5 output format:
            # (0) qName (1) qLength (2) qStart (3) qEnd (4) qStrand
            # (5) tName (6) tLength (7) tStart (8) tEnd (9) tStrand
            # (10...) score ...
            if raw[0] == raw[5]:
                continue  # self-hit
            if raw[4] != raw[9]:
                continue  # opp strand
            f.write(line)

    os.remove(tmp_out)


def pbdagcon_wrapper(fasta_filename, output_prefix,
                     consensus_name, nproc=8,
                     maxScore=-1000, min_seq_len=300):
    """
    (1) Find the best seed as reference
    (2) Align rest to seed
    (3) Call pbdagcon
    """
    ref_filename = output_prefix + '_ref.fasta'
    try:
        out_filename_m1 = output_prefix + ".saln.m1"
        ref = choose_template_by_blasr(fasta_filename=fasta_filename,
                                       out_filename=out_filename_m1,
                                       nproc=nproc, maxScore=maxScore)
        os.remove(out_filename_m1)

        with open(ref_filename, 'w') as f:
            f.write(">{0}\n{1}\n".format(consensus_name, ref.sequence))

        # create alignment file
        aln_filename = output_prefix + '.saln'
        make_aln_input_to_ref(fasta_filename=fasta_filename,
                              ref_filename=ref_filename,
                              out_filename=aln_filename,
                              nproc=nproc)

        cons_filename = output_prefix + '.fasta'
        tmp_cons_filename = output_prefix + '.fasta.tmp'
        # call pbdagcon
        cmd = "pbdagcon -t 0 -m {minlen} -c 1 -j {nproc} {aln} > {out}".format(
            minlen=min_seq_len, nproc=nproc, aln=aln_filename,
            out=tmp_cons_filename)
        cmd += " 2>/dev/null"
        out, code, msg = backticks(cmd)
        if code != 0:
            raise AlignGraphUtilError("Cannot run command: %s" % cmd)

        with FastaReader(tmp_cons_filename) as reader, \
            open(cons_filename, 'w') as writer:
            for rec in reader:
                name = rec.name.strip()
                if "/" in name:
                    # change cid format from c{cid}/0_{len} to c{cid}
                    name = name[:name.find('/')]
                seq = rec.sequence.strip()
                if not 'N' in seq: # Don't write if seq contains N
                    writer.write(">{0}\n{1}\n".format(name, seq))
        os.remove(tmp_cons_filename)

    except AlignGraphUtilError:
        # pick the first sequence as reference as a backup plan
        first_seq = FastaReader(fasta_filename).__iter__().next()
        with open(ref_filename, 'w') as f:
            f.write(">{0}_ref\n{1}\n".
                    format(consensus_name, first_seq.sequence))
    return 0


def set_parser():
    """Set up and return argument parser."""
    parser = ArgumentParser()
    parser.add_argument("input_fasta",
                        help="Input fasta filename")
    parser.add_argument("output_prefix",
                        help="Output filename prefix (ex: g_consensus)")
    parser.add_argument("consensus_id",
                        help="Consensus sequence ID name (ex: consensus)")
    parser.add_argument("--nproc",
                        default=8, type=int,
                        help="Number of processes")
    parser.add_argument("--maxScore", default=-1000, type=int,
                        help="blasr maxScore")
    parser.add_argument("--version", "-v",
                        action='version', version='%(prog)s ' + get_version())
    return parser


def restore_args_with_whitespace(x):
    """Restore args with whitespace, to support file paths with whitespace"""
    y = []
    for xi in x:
        if len(y) == 0:
            y.append(xi)
        else:
            if y[-1].endswith('\\'):
                y[-1] = y[-1][0:-1] + ' ' + xi
            else:
                y.append(xi)
    return y


def runConsensus(args_list_str):
    """Build consesus sequences from fasta input using pbdagcon."""
    parser = set_parser()
    args_list = restore_args_with_whitespace(args_list_str.split())
    args = parser.parse_args(args_list)
    return pbdagcon_wrapper(fasta_filename=args.input_fasta,
                            output_prefix=args.output_prefix,
                            consensus_name=args.consensus_id,
                            nproc=args.nproc, maxScore=args.maxScore)


if __name__ == "__main__":
    sys.exit(runConsensus(" ".join(sys.argv[1:])))