#!/usr/bin/env python
###########################################################################
# Copyright (c) 2011-2014, Pacific Biosciences of California, Inc.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted (subject to the limitations in the
# disclaimer below) provided that the following conditions are met:
#
#  * Redistributions of source code must retain the above copyright
#  notice, this list of conditions and the following disclaimer.
#
#  * Redistributions in binary form must reproduce the above
#  copyright notice, this list of conditions and the following
#  disclaimer in the documentation and/or other materials provided
#  with the distribution.
#
#  * Neither the name of Pacific Biosciences nor the names of its
#  contributors may be used to endorse or promote products derived
#  from this software without specific prior written permission.
#
# NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
# GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
# BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
# USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
# OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
# SUCH DAMAGE.
###########################################################################
"""
Given an input_fasta file of non-full-length (partial) reads and
(unpolished) consensus isoforms sequences in ref_fasta, align reads to
consensus isoforms using BLASR, and then build up a mapping between
consensus isoforms and reads (i.e., assign reads to isoforms).
Finally, save
    {isoform_id: [read_ids],
     nohit: set(no_hit_read_ids)}
to an output pickle file.
"""

import os
import os.path as op
import re
import time
import logging
from pickle import dump
import json

from pbcommand.models import FileTypes
from pbcore.io import ContigSet

from pbtranscript.ClusterOptions import IceOptions
from pbtranscript.Utils import realpath, touch, real_upath, execute
from pbtranscript.PBTranscriptOptions import add_fofn_arguments, \
        add_tmp_dir_argument, add_use_blasr_argument
from pbtranscript.io.ContigSetReaderWrapper import ContigSetReaderWrapper
from pbtranscript.ice_daligner import DalignerRunner

from pbtranscript.ice.IceUtils import set_probqv_from_model, set_probqv_from_fq
from pbtranscript.ice.__init__ import ICE_PARTIAL_PY

from Bio import SeqIO

from cupcake2.ice2.IceUtils2 import blasr_against_ref2, daligner_against_ref2
from cupcake2.tofu2.ToFuOptions2 import add_partial_argument


rex_local_cid = re.compile('c(\d+)')
rex_collected_cid = re.compile('b(\d+)_c(\d+)')

def automatic_determine_if_sID_starts_with_c(seqid_example):
    m = rex_local_cid.match(seqid_example)
    if m is not None:
        return True
    else:
        return False


def build_uc_from_partial_daligner(input_fasta, ref_fasta, out_pickle,
                                   done_filename,
                                   ice_opts,
                                   probqv,
                                   qv_prob_threshold=0.3,
                                   cpus=4,
                                   no_qv_or_aln_checking=False,
                                   tmp_dir=None,
                                   sID_starts_with_c=False):
    """
    Given an input_fasta file of non-full-length (partial) reads and
    (unpolished) consensus isoforms sequences in ref_fasta, align reads to
    consensus isoforms using DALIGNER, and then build up a mapping between
    consensus isoforms and reads (i.e., assign reads to isoforms).
    Finally, save
        {isoform_id: [read_ids],
         nohit: set(no_hit_read_ids)}
    to an output pickle file.

    tmp_dir - where to save intermediate files such as dazz files.
              if None, writer dazz files to the same directory as query/target.
    """
    input_fasta = realpath(input_fasta)
    ref_fasta = realpath(ref_fasta)
    out_pickle = realpath(out_pickle)
    output_dir = op.dirname(out_pickle)

    ice_opts.detect_cDNA_size(ref_fasta)

    # ice_partial is already being called through qsub, so run everything local!
    runner = DalignerRunner(query_filename=input_fasta,
                            target_filename=ref_fasta,
                            is_FL=False, same_strand_only=False,
                            query_converted=False, target_converted=True,
                            dazz_dir=tmp_dir, script_dir=op.join(output_dir, "script"),
                            use_sge=False, sge_opts=None, cpus=cpus)
    runner.run(min_match_len=ice_opts.min_match_len, output_dir=output_dir, sensitive_mode=ice_opts.sensitive_mode)

    partial_uc = {}  # Maps each isoform (cluster) id to a list of reads
    # which can map to the isoform
    seen = set()  # reads seen
    logging.info("Building uc from DALIGNER hits.")

    for la4ice_filename in runner.la4ice_filenames:
        start_t = time.time()


        # not providing full_missed_start/end since aligning nFLs, ok to partially align only
        hitItems = daligner_against_ref2(query_dazz_handler=runner.query_dazz_handler,
                                        target_dazz_handler=runner.target_dazz_handler,
                                        la4ice_filename=la4ice_filename,
                                        is_FL=False, sID_starts_with_c=sID_starts_with_c,
                                        qver_get_func=probqv.get_smoothed,
                                        qvmean_get_func=probqv.get_mean,
                                        qv_prob_threshold=qv_prob_threshold,
                                        ece_penalty=ice_opts.ece_penalty,
                                        ece_min_len=ice_opts.ece_min_len,
                                        same_strand_only=True,
                                        no_qv_or_aln_checking=no_qv_or_aln_checking,
                                        max_missed_start=ice_opts.max_missed_start,
                                        max_missed_end=ice_opts.max_missed_end,
                                        full_missed_start=ice_opts.full_missed_start,
                                        full_missed_end=ice_opts.full_missed_end)


        for h in hitItems:
            if h.ece_arr is not None:
                if h.cID not in partial_uc:
                    partial_uc[h.cID] = set()
                partial_uc[h.cID].add(h.qID)
                seen.add(h.qID)
        logging.info("processing %s took %s sec",
                     la4ice_filename, str(time.time()-start_t))

    for k in partial_uc:
        partial_uc[k] = list(partial_uc[k])

    allhits = set(r.name.split()[0] for r in ContigSetReaderWrapper(input_fasta))

    logging.info("Counting reads with no hit.")
    nohit = allhits.difference(seen)

    logging.info("Dumping uc to a pickle: %s.", out_pickle)
    with open(out_pickle, 'w') as f:
        if out_pickle.endswith(".pickle"):
            dump({'partial_uc': partial_uc, 'nohit': nohit}, f)
        elif out_pickle.endswith(".json"):
            f.write(json.dumps({'partial_uc': partial_uc, 'nohit': nohit}))
        else:
            raise IOError("Unrecognized extension: %s" % out_pickle)

    done_filename = realpath(done_filename) if done_filename is not None \
        else out_pickle + '.DONE'
    logging.debug("Creating %s.", done_filename)
    touch(done_filename)

    # remove all the .las and .las.out filenames
    runner.clean_run()


def _get_fasta_path(file_name):
    if file_name.endswith(".contigset.xml"):
        ds = ContigSet(file_name)
        fasta_files = ds.toExternalFiles()
        assert len(fasta_files) == 1
        return fasta_files[0]
    return file_name


def build_uc_from_partial_blasr(input_fasta, ref_fasta, out_pickle,
                                done_filename,
                                ice_opts,
                                probqv,
                                qv_prob_threshold=0.3,
                                cpus=4,
                                no_qv_or_aln_checking=False,
                                tmp_dir=None,
                                sID_starts_with_c=False):
    """
    Given an input_fasta file of non-full-length (partial) reads and
    (unpolished) consensus isoforms sequences in ref_fasta, align reads to
    consensus isoforms using BLASR, and then build up a mapping between
    consensus isoforms and reads (i.e., assign reads to isoforms).
    Finally, save
        {isoform_id: [read_ids],
         nohit: set(no_hit_read_ids)}
    to an output pickle file.
    """
    input_fasta = _get_fasta_path(realpath(input_fasta))
    m5_file = os.path.basename(input_fasta) + ".blasr"
    if tmp_dir is not None:
        m5_file = op.join(tmp_dir, m5_file)

    out_pickle = realpath(out_pickle)

    cmd = "blasr {i} ".format(i=real_upath(input_fasta)) + \
          "{r} --bestn 100 --nCandidates 200 ".format(r=real_upath(_get_fasta_path(ref_fasta))) + \
          "--nproc {n} -m 5 ".format(n=cpus) + \
          "--maxScore -1000 --minPctIdentity 85 " + \
          "--minAlnLength {a} ".format(a=ice_opts.min_match_len) + \
          "--out {o} ".format(o=real_upath(m5_file)) + \
          "1>/dev/null 2>/dev/null"

    execute(cmd)


    logging.info("Calling blasr_against_ref ...")

    # no need to provide full_missed_start/end for nFLs, since is_FL = False
    hitItems = blasr_against_ref2(output_filename=m5_file,
                                 is_FL=False,
                                 sID_starts_with_c=sID_starts_with_c,
                                 qver_get_func=probqv.get_smoothed,
                                 qvmean_get_func=probqv.get_mean,
                                 qv_prob_threshold=qv_prob_threshold,
                                 ece_penalty=ice_opts.ece_penalty,
                                 ece_min_len=ice_opts.ece_min_len,
                                 max_missed_start=ice_opts.max_missed_start,
                                 max_missed_end=ice_opts.max_missed_end,
                                 full_missed_start=ice_opts.full_missed_start,
                                 full_missed_end=ice_opts.full_missed_end,
                                 same_strand_only=False)


    partial_uc = {}  # Maps each isoform (cluster) id to a list of reads
    # which can map to the isoform
    seen = set()  # reads seen
    logging.info("Building uc from BLASR hits.")
    for h in hitItems:
        if h.ece_arr is not None:
            if h.cID not in partial_uc:
                partial_uc[h.cID] = set()
            partial_uc[h.cID].add(h.qID)
            seen.add(h.qID)

    for k in partial_uc:
        partial_uc[k] = list(partial_uc[k])

    allhits = set(r.name.split()[0] for r in ContigSetReaderWrapper(input_fasta))

    logging.info("Counting reads with no hit.")
    nohit = allhits.difference(seen)

    logging.info("Dumping uc to a pickle: %s.", out_pickle)
    with open(out_pickle, 'w') as f:
        if out_pickle.endswith(".pickle"):
            dump({'partial_uc': partial_uc, 'nohit': nohit}, f)
        elif out_pickle.endswith(".json"):
            f.write(json.dumps({'partial_uc': partial_uc, 'nohit': nohit}))
        else:
            raise IOError("Unrecognized extension: %s" % out_pickle)

    os.remove(m5_file)

    done_filename = realpath(done_filename) if done_filename is not None \
        else out_pickle + '.DONE'
    logging.debug("Creating %s.", done_filename)
    touch(done_filename)


class IcePartialOne2(object):

    """Assign nfl reads of a given fasta to isoforms."""

    desc = "Assign non-full-length reads in the given input fasta to " + \
           "unpolished consensus isoforms."
    prog = "%s one " % ICE_PARTIAL_PY

    def __init__(self, input_fasta, input_fastq, ref_fasta, out_pickle,
                 done_filename, ice_opts, cpus=4, tmp_dir=None):
        self.input_fasta = input_fasta
        self.input_fastq = input_fastq # could be None
        self.ref_fasta = ref_fasta
        self.out_pickle = out_pickle
        self.done_filename = done_filename
        self.cpus = cpus
        self.ice_opts = ice_opts
        self.tmp_dir = tmp_dir

        # read QV from fastq
        start_t = time.time()
        if self.input_fastq is not None:
            self.probqv, msg = set_probqv_from_fq(self.input_fastq)
        else:
            self.probqv, msg = set_probqv_from_model()
        logging.info("Reading probQV from {0} took {1:.1f} sec.".format(self.input_fastq, time.time()-start_t))


    def cmd_str(self):
        """Return a cmd string (ice_partial.py one)."""
        return self._cmd_str(input_fasta=self.input_fasta,
                             ref_fasta=self.ref_fasta,
                             out_pickle=self.out_pickle,
                             done_filename=self.done_filename,
                             ice_opts=self.ice_opts,
                             cpus=self.cpus,
                             tmp_dir=self.tmp_dir)

    def _cmd_str(self, input_fasta, ref_fasta, out_pickle,
                 done_filename, ice_opts, cpus, tmp_dir=None):
        """Return a cmd string (ice_partil.py one)"""
        cmd = self.prog + \
              "{f} ".format(f=input_fasta) + \
              "{r} ".format(r=ref_fasta) + \
              "{o} ".format(o=out_pickle)
        if done_filename is not None:
            cmd += "--done {d} ".format(d=done_filename)
        cmd += ice_opts.cmd_str()
        cmd += "--cpus={0}".format(cpus)
        if tmp_dir is not None:
            cmd += "--tmp_dir {t} ".format(t=tmp_dir)
        return cmd



    def run(self):
        """Run"""
        # automatically determine whether the sIDs are local
        # ex: c0/123
        # or after collection
        # ex: b<bin>_c<cid>

        r = next(SeqIO.parse(open(self.ref_fasta),'fasta'))
        sID_starts_with_c = automatic_determine_if_sID_starts_with_c(r.id)
        logging.info("sID_starts_with_c: {0}.".format(sID_starts_with_c))

        logging.info("Building uc from non-full-length reads using {0}.".format(self.ice_opts.aligner_choice))



        if self.ice_opts.aligner_choice == 'daligner':
            try:
                build_uc_from_partial_daligner(input_fasta=self.input_fasta,
                                               ref_fasta=self.ref_fasta,
                                               out_pickle=self.out_pickle,
                                               done_filename=self.done_filename,
                                               ice_opts=self.ice_opts,
                                               probqv=self.probqv,
                                               qv_prob_threshold=0.3,
                                               cpus=self.cpus,
                                               no_qv_or_aln_checking=False,
                                               tmp_dir=self.tmp_dir,
                                               sID_starts_with_c=sID_starts_with_c)
            except:
                # fall back to blasr
                self.ice_opts.aligner_choice = 'blasr'
                self.run()
        elif self.ice_opts.aligner_choice == 'blasr':
            build_uc_from_partial_blasr(input_fasta=self.input_fasta,
                                  ref_fasta=self.ref_fasta,
                                  out_pickle=self.out_pickle,
                                  done_filename=self.done_filename,
                                  ice_opts=self.ice_opts,
                                  probqv=self.probqv,
                                  qv_prob_threshold=0.3,
                                  cpus=self.cpus,
                                  no_qv_or_aln_checking=False,
                                  tmp_dir=self.tmp_dir,
                                  sID_starts_with_c=sID_starts_with_c)
        else:
            raise Exception("Aligner choice {0} not recognized!".format(self.ice_opts.aligner_choice))
        return 0


def add_ice_partial_one_arguments(parser):
    """Add arguments for assigning nfl reads of a given input fasta
    to unpolished isoforms."""
    parser.add_argument("input_fasta", help="Non full-length reads (fasta)")
    parser.add_argument("--input_fastq", default=None, help="Non full-length reads (fastq). If not given, no QVs are used.")
    parser.add_argument("ref_fasta", help="Reference fasta, most likely ref.consensus.fasta")
    parser.add_argument("out_pickle", help="Output pickle file", default="ice_partial_one")
    parser.add_argument("--done", dest="done_filename", type=str,
                            help="An empty file generated to indicate that " +
                            "out_pickle is done.")
    parser = add_partial_argument(parser)
    parser = add_tmp_dir_argument(parser)


    return parser