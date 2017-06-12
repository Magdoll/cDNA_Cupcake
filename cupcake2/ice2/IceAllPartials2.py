#!/usr/bin/env python
"""
IceAllPartials procedures:

(1) For each input fasta file (i.e., each splitted nfl read file), map its
    reads to unpolished consensus isoforms in ref_fasta
    (i.e., final.consensus.fasta), then create a pickle file
    (e.g. *.partial_uc.pickle).

(2) Wait for all pickle files to be created

(3) Merge all pickle files and dump to a big pickle.

"""

import os.path as op
import logging
import time
#from pbtranscript.PBTranscriptOptions import \
#    add_sge_arguments, add_fofn_arguments, add_tmp_dir_argument
from pbtranscript.Utils import realpath, mkdir, real_upath, ln
from pbtranscript.ice.IceFiles import IceFiles
from pbtranscript.ice.IceUtils import combine_nfl_pickles


from cupcake2.ice2.__init__ import ICE_PARTIAL_PY


from cupcake2.tofu2.ToFuOptions2 import add_sge_arguments, add_fofn_arguments, add_tmp_dir_argument, add_partial_argument

class IceAllPartials2(IceFiles):

    """Align non-full-length reads to unpolished consensus
    isoforms and assign them to isoforms based on blasr alignments,
    (build partial uc) and save results to a pickle file.
    """

    # Define description of IceAllPartials
    desc = "Align and assign all non-full-length reads to unpolished " + \
           "consensus isoforms, and save results to a pickle file."

    prog = "%s all " % ICE_PARTIAL_PY


    def __init__(self, root_dir, fasta_filenames, fastq_filenames, ref_fasta,
                 out_pickle, ice_opts, sge_opts, cpus, tmp_dir=None):
        """
        fasta_filenames --- a list of splitted nfl fasta files.

        ref_fasta --- (unpolished) consensus isoforms

        out_pickle --- a pickle file with all nfl fasta reads

        root_dir --- ICE root output directory

        tmp_dir --- if not None, write temporary clusters, dazz, las
                    files to the given temporaray directory

        sge_opts --- params for SGE environment, including
            use_sge    : use SGE or not
            max_sge_jobs: maximum number of sub-jobs submitted
            unique_id  : unique qsub job id, important that this
                        DOES NOT CONFLICT!

        """
        self.prog_name = "IceAllPartials"
        IceFiles.__init__(self, prog_name=self.prog_name, root_dir=root_dir, tmp_dir=tmp_dir)

        self.fasta_filenames, self.ref_fasta = \
            self._validate_inputs(fasta_filenames=fasta_filenames,
                                  ref_fasta=ref_fasta)
        self.fastq_filenames = fastq_filenames

        self.out_pickle = out_pickle

        self.ice_opts = ice_opts
        self.sge_opts = sge_opts
        self.cpus = cpus # this is the number of CPUs to use per SGE job or per local job


        self.add_log("Making dir for mapping noFL reads: " + self.nfl_dir)
        mkdir(self.nfl_dir)

        self.add_log("input fasta files are: " +
                     ", ".join(self.fasta_filenames))
        self.add_log("temp pickle files are: " +
                     ", ".join(self.pickle_filenames))
        self.add_log("out pickle file is: " + self.out_pickle)
        self.add_log("temp directory is: " + str(self.tmp_dir))

    def _validate_inputs(self, fasta_filenames, ref_fasta):
        """Validate input files."""
        for f in fasta_filenames:
            if not op.exists(f):
                raise IOError("Input fasta {f} does not exist.".format(f=f))
        if ref_fasta is None or not op.exists(ref_fasta):
            raise IOError("Reference {r} does not exist.".format(r=ref_fasta))

        return ([realpath(f) for f in fasta_filenames],
                realpath(ref_fasta))

    def cmd_str(self):
        """Return a cmd string."""
        return self._cmd_str(fasta_filenames=self.fasta_filenames,
                             ref_fasta=self.ref_fasta,
                             out_pickle=self.out_pickle,
                             root_dir=self.root_dir,
                             ice_opts=self.ice_opts,
                             sge_opts=self.sge_opts,
                             cpus=self.cpus,
                             tmp_dir=self.tmp_dir)

    def _cmd_str(self, fasta_filenames, ref_fasta, out_pickle,
                 root_dir, ice_opts, sge_opts, cpus, tmp_dir):
        """Return a cmd string."""
        cmd = self.prog + \
              "{fns} ".format(fns=",".join(fasta_filenames)) + \
              "{ref} ".format(ref=ref_fasta) + \
              "{out} ".format(out=out_pickle)
        if root_dir is not None:
            cmd += "--root_dir={d} ".format(d=root_dir)
        cmd += "--cpus={0}".format(cpus)
        cmd += ice_opts.cmd_str()
        cmd += sge_opts.cmd_str()
        if tmp_dir is not None:
            cmd += "--tmp_dir={t} ".format(t=tmp_dir)
        return cmd

    @property
    def pickle_filenames(self):
        """pickle files for each fasta file."""
        return [op.join(self.nfl_dir, op.basename(f) + ".partial_uc.pickle")
                for f in self.fasta_filenames]

    @property
    def done_filenames(self):
        """done files to indicate that pickles are done."""
        return [op.join(self.nfl_dir, op.basename(f) + ".DONE")
                for f in self.pickle_filenames]

    @property
    def script_filenames(self):
        """scripts to generate pickles from fasta files."""
        return [op.join(self.script_dir, op.basename(f) + ".partial_uc.sh")
                for f in self.fasta_filenames]

    def createPickles(self):
        """For each file in fasta_filenames, call 'ICE_PARTIAL_PY one' to
        build clusters and to save results to a pickle file. When all pickles
        are done, union all pickles.
        """
        self.add_log("Mapping non-full-length reads to consensus isoforms.")
        self.add_log("Creating pickles...", level=logging.INFO)

        for idx, fa in enumerate(self.fasta_filenames):
            # for each splitted non-full-length reads fasta file, build #
            # partial_uc.pickle

            # ex:
            #      python run_IcePartial2.py one isoseq_nfl.fasta isoseq_nfl.fastq \
            #  output/final.consensus.fasta isoseq_nfl.fasta.pickle --aligner_choice=blasr  --cpus=12
            fq = self.fastq_filenames[idx]

            cmd = ICE_PARTIAL_PY + " " + \
                  "one {fa} {fq} ".format(fa=real_upath(fa), fq=real_upath(fq)) + \
                  "{r} ".format(r=real_upath(self.ref_fasta)) + \
                  "{o} ".format(o=real_upath(self.pickle_filenames[idx])) + \
                  "--aligner_choice={c} ".format(c=self.ice_opts.aligner_choice) + \
                  "--cpus={n} ".format(n=self.sge_opts.blasr_nproc) + \
                  "--max_missed_start={0} ".format(self.ice_opts.max_missed_start) + \
                  "--max_missed_end={0} ".format(self.ice_opts.max_missed_end) + \
                  "--ece_penalty={0} ".format(self.ice_opts.ece_penalty) + \
                  "--ece_min_len={0} ".format(self.ice_opts.ece_min_len) + \
                  "--done={d} ".format(d=real_upath(self.done_filenames[idx]))
            if self.tmp_dir is not None:
                cmd += "--tmp_dir={t}".format(t=self.tmp_dir)

            self.add_log("Writing command to script {fsh}".
                         format(fsh=self.script_filenames[idx]))
            with open(self.script_filenames[idx], 'w') as fsh:
                fsh.write(cmd + "\n")

            # determine elog & olog
            partial_log_fn = op.join(self.log_dir,
                                     'IcePartial.{idx}'.format(idx=idx))
            elog = partial_log_fn + ".elog"
            olog = partial_log_fn + ".olog"
            jid = "ice_partial_{unique_id}_{name}".format(
                unique_id=self.sge_opts.unique_id,
                name=op.basename(fa))

            self.add_log("Creating a pickle for {f}".format(f=fa))

            if self.sge_opts.use_sge is True:
                qsub_cmd = self.sge_opts.qsub_cmd(script=real_upath(self.script_filenames[idx]),
                                                  num_threads=self.cpus,
                                                  wait_before_exit=False,
                                                  depend_on_jobs=None,
                                                  elog=real_upath(elog),
                                                  olog=real_upath(olog),
                                                  is_script=True,
                                                  jobid=jid)
                #          qsub_cmd = "qsub " + \
                #                     "-pe smp {n} ".format(n=self.sge_opts.blasr_nproc) + \
                #                     "-cwd -S /bin/bash -V " + \
                #                     "-e {elog} ".format(elog=real_upath(elog)) + \
                #                     "-o {olog} ".format(olog=real_upath(olog)) + \
                #                     "-N {jid} ".format(jid=jid) + \
                #                     "{sh}".format(sh=real_upath(self.script_filenames[idx]))
                self.qsub_cmd_and_log(qsub_cmd)
            else:
                cmd += " 1>{olog} 2>{elog}".format(olog=real_upath(olog),
                                                   elog=real_upath(elog))
                self.run_cmd_and_log(cmd=cmd, olog=olog, elog=elog)

    def waitForPickles(self, pickle_filenames, done_filenames):
        """Wait for *.pickle and *.pickle.DONE to be created."""
        self.add_log("Waiting for pickles {ps} to be created.".
                     format(ps=", ".join(pickle_filenames)),
                     level=logging.INFO)
        stop = False
        sleep_time = 10
        while stop is not True:
            stop = all(op.exists(p) for p in pickle_filenames) and \
                all(op.exists(d) for d in done_filenames)
            sleep_time = min(600, sleep_time * 2)
            time.sleep(sleep_time)
            self.add_log("Waiting for pickles to be created: {ps}".
                         format(ps=", ".join([p for p in pickle_filenames
                                              if op.exists(p)])))

    def combinePickles(self, pickle_filenames, out_pickle):
        """Combine all *.pickle files to one and dump to self.out_pickle."""
        combine_nfl_pickles(pickle_filenames, out_pickle)

    def run(self):
        """Assigning nfl reads to consensus isoforms and merge."""
        # Call $ICE_PARTIAL_PY to create a pickle for each splitted nfl fasta
        self.createPickles()
        # Wait for pickles to be created, if SGE is used.
        self.waitForPickles(pickle_filenames=self.pickle_filenames,
                            done_filenames=self.done_filenames)
        # Combine all pickles to a big pickle file: nfl_all_pickle_fn.
        self.combinePickles(pickle_filenames=self.pickle_filenames,
                            out_pickle=self.nfl_all_pickle_fn)
        # Create symbolic link if necessary
        ln(self.nfl_all_pickle_fn, self.out_pickle)

        # Close log
        self.close_log()



def add_ice_all_partials_arguments(parser):
    """Add IceAllPartials argument parser."""
    parser.add_argument("fasta_filenames",
                        type=str,
                        help="comma delimited fasta files of " +
                             "splitted non-full-length reads")
    parser.add_argument("ref_fasta",
                        type=str,
                        help="Reference fasta file, most likely " +
                             "final.consensus.fasta from ICE output")
    parser.add_argument("out_pickle",
                        type=str,
                        help="Output pickle file.")
    parser.add_argument("--fastq_filenames",
                        default='',
                        type=str,
                        help="comma delimited fastq files of " +
                             "splitted non-full-length reads")
    parser.add_argument("--root_dir",
                        dest="root_dir",
                        default="",
                        help="A directory for saving intermediate files.")

    parser = add_partial_argument(parser)
    # for running Partial, only relevant nproc is the gcon_nproc
    parser = add_sge_arguments(parser, blasr_nproc=False, arrow_nproc=False, gcon_nproc=False)
    # Liz: I don't think fofn (subread XML) is needed at partial step. commenting out for now.
    #parser = add_fofn_arguments(parser, bas_fofn=True, fasta_fofn=False)

    parser.add_argument("--done", dest="done_filename", type=str,
                        help="An empty file generated to indicate that " +
                             "out_pickle is done.")

    parser = add_tmp_dir_argument(parser)
    return parser
