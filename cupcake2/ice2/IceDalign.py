#!/usr/bin/env python
"""
Class DalignerRunner which aligns a query FASTA file and
a target FASTA file using daligner and interpret las output
using LA4Ice.

daligner jobs can either be submitted to SGE or run locally.

"""

import os
import os.path as op
import sys
import logging
import time

from pbcore.util.ToolRunner import PBToolRunner

from pbtranscript.__init__ import get_version
from pbtranscript.PBTranscriptOptions import add_sge_arguments
from pbtranscript.ClusterOptions import SgeOptions
from pbtranscript.Utils import realpath, mkdir, mknewdir
from pbtranscript.RunnerUtils import write_cmd_to_script, \
    sge_job_runner, local_job_runner
from pbtranscript.io import DazzIDHandler

__author__ = "etseng@pacificbiosciences.com"

#logger = logging.getLogger(op.basename(__file__))
#logger.setLevel(logging.DEBUG)

#NTHREADS is hard-coded as 4 in daligner.
DALIGNER_NUM_THREADS = 4

class DalignerRunner(object):
    """
    DalignerRunner, which aligns query FASTA file to target
    FASTA files using daligner and interpret las output using LA4ice.
    """
    def __init__(self, query_filename, target_filename,
                 is_FL, same_strand_only,
                 query_converted=False, target_converted=False,
                 dazz_dir=None, script_dir="scripts/",
                 use_sge=False, sge_opts=None, cpus=24):
        """
        Parameters:
          query_filename - query FASTA file
          target_filename - target FASTA file
          is_FL - whether or not reads are FLNC CCS reads
          same_strand_only - whether or not align reads in reverse strand
          query_converted - whether or not query FASTA file been converted
                            to daligner compatible FASTA file.
          target_converted - whether or not target FASTA file been converted
                             to daligner compatible FASTA file.
          dazz_dir - if None, all query.dazz.* files will be saved in the same
                     directory as query and all target.dazz.* files will be
                     saved in the same dir as target.
                     if a valid path, all query.dazz.* files and target.dazz.*
                     files will be saved to dazz_dir.
          script_dir - directory for saving all scripts

          use_sge - submit daligner jobs to sge or run them locally?
          sge_opts - sge options
          cpus - total number of cpus that can be used to align query to target.
        """
        self.query_filename = realpath(query_filename)
        self.target_filename = realpath(target_filename)
        self.is_FL = is_FL
        self.same_strand_only = same_strand_only
        self.cpus = cpus
        self.dazz_dir = dazz_dir
        self.script_dir = realpath(script_dir)
        self.output_dir = ""

        self.query_dazz_handler = DazzIDHandler(self.query_filename,
                                                converted=query_converted,
                                                dazz_dir=dazz_dir)
        # target may have already been converted (if shared)
        target_converted = (target_converted or
                            self.query_filename == self.target_filename)
        self.target_dazz_handler = DazzIDHandler(self.target_filename,
                                                 converted=target_converted,
                                                 dazz_dir=dazz_dir)

        self.target_blocks = self.target_dazz_handler.num_blocks
        self.query_blocks = self.query_dazz_handler.num_blocks

        self.use_sge = use_sge
        self.sge_opts = sge_opts

    def query_prefix(self, i):
        """Return (possibly absolute path) prefix of query block i.
        e.g., query.dazz.fasta.{i} if query_blocks > 1
              query.dazz.fasta if query_blocks == 1
        """
        assert i >= 1 and i <= self.query_blocks
        return self.query_dazz_handler.dazz_filename + \
                (".{i}".format(i=i) if self.query_blocks != 1 else "")

    def target_prefix(self, j):
        """Return (possibly absolute path) prefix of target block j.
        e.g., target.dazz.fasta.{j} if target_blocks > 1
              target.dazz.fasta if target_blocks == 1
        """
        assert j >= 1 and j <= self.target_blocks
        return self.target_dazz_handler.dazz_filename + \
                (".{j}".format(j=j) if self.target_blocks != 1 else "")

    def thread_prefix(self, k, is_forward):
        """Return prefix of thread k.
        e.g., N{k} if is forward strand.
        e.g., C{k} if is backward strand.
        """
        assert k >= 0 and k <= DALIGNER_NUM_THREADS - 1
        strand = "N" if is_forward else "C"
        return "{strand}{k}".format(strand=strand, k=k)

    def _iter_i_j(self):
        """Iterate over indices (query block i, target block j)."""
        for i in xrange(1, self.query_blocks + 1):
            for j in xrange(1, self.target_blocks + 1):
                if not (self.query_filename == self.target_filename and i > j):
                    yield (i, j)

    def _iter_i_j_k(self):
        """Iterate over indices (query block i, target block j, thread k)."""
        for k in xrange(DALIGNER_NUM_THREADS):
            for i, j in self._iter_i_j():
                yield (i, j, k)

    def las_filename(self, i, j, k, is_forward, output_dir, switch_query_target=False):
        """Return las output for query block i, target block j, thread k.
        e.g.,query.dazz.fasta.{i}.basename(target).dazz.fasta.{j}.N{k}.las
        """
        _query = op.basename(self.query_prefix(i=i))
        _target = op.basename(self.target_prefix(j=j))
        if switch_query_target is True:
            _target = op.basename(self.query_prefix(i=i))
            _query = op.basename(self.target_prefix(j=j))
        return op.join(output_dir,
                       "{query_prefix}.{target_prefix}.{thread_prefix}.las"
                       .format(query_prefix=_query, target_prefix=_target,
                               thread_prefix=self.thread_prefix(k=k, is_forward=is_forward)))

    def _las_filenames(self, output_dir, is_forward, switch_query_target=False):
        """Return las output file name in a list."""
        return [self.las_filename(i=i, j=j, k=k, is_forward=is_forward,
                                  output_dir=output_dir,
                                  switch_query_target=switch_query_target)
                for i, j, k in self._iter_i_j_k()]

    @property
    def las_filenames(self):
        """Return all las output file name in a list."""
        ret = self._las_filenames(is_forward=True, output_dir=self.output_dir)
        if not self.same_strand_only:
            ret.extend(self._las_filenames(is_forward=False, output_dir=self.output_dir))
        return ret

    def la4ice_filename(self, i, j, k, is_forward, output_dir):
        """Return la4ice output for query block i, target block j, thread k.
        e.g., op.join(output_dir, {las_filenam}.out)
        """
        return op.join(output_dir,
                       self.las_filename(i, j, k, is_forward=is_forward,
                                         output_dir=output_dir) + ".out")

    def _la4ice_filenames(self, is_forward, output_dir):
        """Return a list of la4ice output files showing either forward-only
        or reverse-only alignments"""
        return [self.la4ice_filename(i=i, j=j, k=k, is_forward=is_forward, output_dir=output_dir)
                for i, j, k in self._iter_i_j_k()]

    @property
    def la4ice_filenames(self):
        """Return all la4ice output files as a list of strings."""
        ret = self._la4ice_filenames(is_forward=True, output_dir=self.output_dir)
        if not self.same_strand_only:
            ret.extend(self._la4ice_filenames(is_forward=False, output_dir=self.output_dir))
        return ret

    def daligner_cmd(self, i, j,
                     min_match_len=300, sensitive_mode=False):
        """Return dalianger command running for query block i, target block j.
        e.g., daligner -h35 -e.80 -l500 -s100 <query> <target>
        Note that daligner always writes output to the current dir, so we
        need to change workding directory when running the command later.
        """
        # old DALIGNER param is not sensitive enough for > 5 kb CCS reads
        # but the more sensitive param seems to do badly on 1 - 2 kb reads. WTH =_=
        prog = "daligner"
        params = "-h35 -k16 -e.80 -l{m} -s100 -t10 ".format(m=min_match_len)
        if sensitive_mode:
            params = "-w12 -h24 -k24 -e.70 -l{m} -s100 -t10 ".format(m=min_match_len)
        return " ".join([prog, params, self.query_prefix(i=i), self.target_prefix(j=j)])

    def daligner_cmds(self, min_match_len=300, sensitive_mode=False):
        """Return all daligner commands as a list."""
        # avoid unnecessary redundant calls to daligner when doing self hits
        return [self.daligner_cmd(i=i, j=j, min_match_len=min_match_len,
                                  sensitive_mode=sensitive_mode)
                for i, j in self._iter_i_j()]

    @property
    def daligner_scripts(self):
        """Return daligner scripts."""
        return [op.join(self.script_dir, "daligner_{i}_{j}.sh".format(i=i, j=j))
                for i, j in self._iter_i_j()]

    def la4ice_cmd(self, i, j, k, is_forward, output_dir):
        """Return a la4ice command for query block i, target block j."""
        prog = "LA4Ice"
        params = "-a -m -i0 -w100000 -b0 "
        if self.is_FL:
            params += "-E"
        return " ".join([prog, params, self.query_dazz_handler.dazz_filename,
                         self.target_dazz_handler.dazz_filename,
                         self.las_filename(i=i, j=j, k=k, is_forward=is_forward,
                                           output_dir=output_dir),
                         ">",
                         self.la4ice_filename(i=i, j=j, k=k, is_forward=is_forward,
                                              output_dir=output_dir)])

    def _la4ice_cmds(self, is_forward, output_dir):
        """Return a list of la4ice commands showing either forward only or
        reverse only alignments.
        """
        return [self.la4ice_cmd(i=i, j=j, k=k, is_forward=is_forward, output_dir=output_dir)
                for i, j, k in self._iter_i_j_k()]

    @property
    def la4ice_cmds(self):
        """Return la4ice commands as a list of strings"""
        ret = self._la4ice_cmds(is_forward=True, output_dir=self.output_dir)
        if not self.same_strand_only:
            ret.extend(self._la4ice_cmds(is_forward=False, output_dir=self.output_dir))
        return ret

    @property
    def la4ice_scripts(self):
        """Return la4ice scripts."""
        ret = [op.join(self.script_dir, "la4ice_{i}_{j}_N{k}.sh".format(i=i, j=j, k=k))
               for i, j, k in self._iter_i_j_k()]
        if not self.same_strand_only:
            ret.extend([op.join(self.script_dir, "la4ice_{i}_{j}_C{k}.sh".format(i=i, j=j, k=k))
                        for i, j, k in self._iter_i_j_k()])
        return ret

    @property
    def daligner_done_script(self):
        """Return script_dir/daligner_done.sh"""
        return op.join(self.script_dir, "daligner_done.sh")

    def write_daligner_done_script(self, output_dir):
        """
        Write cmd to self.daligner_done_script, the cmd creates
        a sentinel file indicating all daligner jobs completed.

        The sentinel file that will be written is output_dir/DALIGNER.DONE.

        """
        cmd = "touch {out}".format(out=op.join(output_dir, "DALIGNER.DONE"))
        write_cmd_to_script(cmd=cmd, script=self.daligner_done_script)

    @property
    def la4ice_done_script(self):
        """Return script_dir/la4ice_done.sh"""
        return op.join(self.script_dir, "la4ice_done.sh")

    def write_la4ice_done_script(self, output_dir):
        """
        Write cmd to self.la4ice_done_script, the cmd creates
        a sentinel file indicating all la4ice jobs compeleted.

        The sentinel file that will be written is output/LA4ICE.DONE.
        """
        cmd = "touch {out}".format(out=op.join(output_dir, "LA4ICE.DONE"))
        write_cmd_to_script(cmd=cmd, script=self.la4ice_done_script)

    def la4ice_script(self, i):
        """Return script_dir/la4ice_job_{i}.sh"""
        return op.join(self.script_dir, "LA4Ice_job_{i}.sh".format(i=i))

    def run(self, output_dir='.', min_match_len=300, sensitive_mode=False):
        """
        if self.use_sge --- writes to <scripts>/daligner_job_#.sh
        else --- run locally, dividing into self.cpus/4 tasks (capped max at 4)

        NOTE 1: when using SGE, be careful that multiple calls to this might
        end up writing to the SAME job.sh files, this should be avoided by
        changing <scripts> directory

        NOTE 2: more commonly this should be invoked locally
        (since ice_partial.py i/one be qsub-ed),
        in that case it is more recommended to keep self.cpus = 4 so that
        each daligner job is run consecutively and that the original qsub job
        should have been called with qsub -pe smp 4 (set by --blasr_nproc 4)
        In this way, the daligner jobs are called consecutively, but LA4Ice
        is parallelized 4X
        """
        self.output_dir = realpath(output_dir) # Reset output_dir
        old_dir = realpath(op.curdir)
        mkdir(output_dir)
        os.chdir(output_dir)

        if self.use_sge:
            mknewdir(self.script_dir)

        # prepare done scripts is no longer necessary.
        #self.write_daligner_done_script()
        #self.write_la4ice_done_script()

        # (a) run all daligner jobs
        daligner_cmds = self.daligner_cmds(min_match_len=min_match_len,
                                           sensitive_mode=sensitive_mode)

        logging.info("Start daligner cmds " +
                     ("using sge." if self.use_sge else "locally."))
        logging.debug("CMD: " + "\n".join(daligner_cmds))

        start_t = time.time()
        failed = []
        if self.use_sge:
            failed.extend(
                sge_job_runner(cmds_list=daligner_cmds,
                               script_files=self.daligner_scripts,
                               #done_script=self.daligner_done_script,
                               num_threads_per_job=DALIGNER_NUM_THREADS,
                               sge_opts=self.sge_opts, qsub_try_times=3,
                               wait_timeout=600, run_timeout=600,
                               rescue="sge", rescue_times=3))
        else:
            # max 4 at a time to avoid running out of memory...
            failed.extend(
                local_job_runner(cmds_list=daligner_cmds,
                                 num_threads=max(1, min(self.cpus/4, 4))))
        logging.info("daligner jobs took " + str(time.time()-start_t) + " sec.")

        # (b) run all LA4Ice jobs
        start_t = time.time()
        logging.info("Start LA4Ice cmds " +
                     ("using sge." if self.use_sge else "locally."))
        la4ice_cmds = self.la4ice_cmds
        logging.debug("CMD: " + "\n".join(la4ice_cmds))

        if self.use_sge:
            failed.extend(
                sge_job_runner(cmds_list=la4ice_cmds,
                               script_files=self.la4ice_scripts,
                               #done_script=self.la4ice_done_script,
                               num_threads_per_job=DALIGNER_NUM_THREADS,
                               sge_opts=self.sge_opts, qsub_try_times=3,
                               wait_timeout=600, run_timeout=600,
                               rescue="sge", rescue_times=3))
        else:
            # max 4 at a time to avoid running out of memory...
            failed.extend(
                local_job_runner(cmds_list=la4ice_cmds,
                                 num_threads=max(1, min(self.cpus, 4))))
        logging.info("LA4Ice jobs took " + str(time.time()-start_t) + " sec.")
        os.chdir(old_dir)

        if len(failed) == 0:
            return 0
        else:
            raise RuntimeError("%s.run failed, %s." %
                               (op.basename(self.__class__),
                                "\n".join([x[0] for x in failed])))

    def clean_run(self):
        """Remove output files: las_filenames, la4ice_filenames."""
        fs = self._las_filenames(is_forward=True, output_dir=self.output_dir, switch_query_target=False) + \
             self._las_filenames(is_forward=False, output_dir=self.output_dir, switch_query_target=False) + \
             self._las_filenames(is_forward=True, output_dir=self.output_dir, switch_query_target=True) + \
             self._las_filenames(is_forward=False, output_dir=self.output_dir, switch_query_target=True) + \
             self.la4ice_filenames
        for f in set(fs):
            os.remove(f)


def add_ice_daligner_arguments(parser):
    """Set up argument parser."""
    helpstr = "Query reads fasta file"
    parser.add_argument("query_fasta", type=str, help=helpstr)

    helpstr = "Target reads fasta file"
    parser.add_argument("target_fasta", type=str, help=helpstr)

    helpstr = "Output directory to store daligner and LA4Ice outputs"
    parser.add_argument("output_dir", type=str, help=helpstr)

    helpstr = "Query reads are Full-Length isoforms."
    parser.add_argument("--is_FL", default=False, action="store_true", help=helpstr)

    helpstr = "Query and target reads are of the same strand"
    parser.add_argument("--same_strand_only", default=False, action="store_true", help=helpstr)

    parser = add_sge_arguments(parser, blasr_nproc=True)
    return parser


class IceDalignerRunner(PBToolRunner):

    "Use daligner to align cDNA reads, then call LA4Ice to show alignments."

    def __init__(self):
        desc = __doc__
        PBToolRunner.__init__(self, desc)
        add_ice_daligner_arguments(self.parser)

    def getVersion(self):
        """Get version string."""
        return get_version()

    def run(self):
        """ Call DalignerRunner """
        logging.info("Running {f} v{v}.".format(f=op.basename(__file__),
                                                v=self.getVersion()))
        args = self.args
        mkdir(args.output_dir)

        sge_opts = SgeOptions(unique_id=args.unique_id,
                              use_sge=args.use_sge,
                              max_sge_jobs=args.max_sge_jobs,
                              blasr_nproc=args.blasr_nproc,
                              sge_env_name=args.sge_env_name,
                              sge_queue=args.sge_queue)

        obj = DalignerRunner(query_filename=args.query_fasta,
                             target_filename=args.target_fasta,
                             is_FL=args.is_FL, same_strand_only=args.same_strand_only,
                             query_converted=False, target_converted=False,
                             use_sge=args.use_sge, sge_opts=sge_opts)
        obj.run(output_dir=args.output_dir)


def main():
    """Main function to call CombineClusterBinsRunner"""
    runner = IceDalignerRunner()
    return runner.start()


if __name__ == "__main__":
    sys.exit(main())