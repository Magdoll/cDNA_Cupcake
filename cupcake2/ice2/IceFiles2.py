__author__ = 'etseng@pacb.com'


"""
Diff with IceFiles:
  1. <subread_xml> replaces <bas_fofn> (for running Arrow)
  2. removed <fasta_fofn>
  3. removed <ccs_fofn>
  4. s/quiver/arrow
"""

import os.path as op
import logging
from pbtranscript.Utils import real_ppath, now_str, mkdir

from pbtranscript.ice.IceFiles import IceFiles

class IceFiles2(IceFiles):
    nfl_fa_format = "input.split_{0:03d}.fasta"
    nfl_fq_format = "input.split_{0:03d}.fastq"


    def __init__(self, prog_name, root_dir,
                 subread_xml=None,
                 no_log_f=False,
                 tmp_dir=None,
                 make_dirs=True):
        """
        prog_name --- name of a sub-class
        root_dir --- root directory of the whole project. There will be
                     sub-directories under it, including:
                     tmp/ --- 0/  c0, c1, ..., c9999
                          --- 1/  c10000, c10001, ..., c19999
                          ...
                          each c? folder contains data for a cluster id=c?
                     script/
                          --- 0/  gcon_job_?.sh, gcon jobs in the first iteration
                          --- 1/  gcon_job_?.sh, gcon jobs in the second iteration
                          ...
                     log/
                          --- ICE.log   Log of the ICE algorithm
                          --- 0/  log for jobs in the first iteration
                          ...
                     output/   output files go here.
        subread_xml --- subreads.xml file
        no_log_f --- DON'T write log to a log file.
        tmp_dir --- Write temporary files to tmp_dir (usually /scratch) for speed
        """
        self.prog_name = str(prog_name)
        self.root_dir = real_ppath(root_dir)
        self._tmp_dir = real_ppath(tmp_dir)

        self.subread_xml = real_ppath(subread_xml)

        if make_dirs is True:
            mkdir(self.root_dir)
            mkdir(self.tmp_dir)
            mkdir(self.log_dir)
            mkdir(self.script_dir)
            mkdir(self.out_dir)

        self.no_log_f = no_log_f
        if not no_log_f:
            self.log_f = open(self.log_fn, 'w', 0)
            self.add_log(msg="{p} initialized.".format(p=self.prog_name))


    def nfl_fa_i(self, i):
        return op.join(self.nfl_dir, IceFiles2.nfl_fa_format.format(i))

    def nfl_fq_i(self, i):
        """Return the i-th splitted chunk of nfl reads.
           $root_dir/output/map_noFL/input.split_{0:03d}.fastq
        """
        return op.join(self.nfl_dir, IceFiles2.nfl_fq_format.format(i))


    @property
    def arrowed_dir(self):
        """Return $root_dir/arrowed"""
        return op.join(self.root_dir, "arrowed")

    @property
    def arrowed_log_dir(self):
        """Return $log_dir/arrowed"""
        return op.join(self.log_dir, "arrowed")

    @property
    def prepare_arrow_files(self):
        """Return $root_dir/log/prepared_arrow_files.txt"""
        return op.join(self.log_dir, 'prepared_arrow_files.txt')

    def arrow_submission_file(self, i, total):
        return op.join(self.log_dir, "arrow_job_{0}of{1}.sh".format(i, total))

    @property
    def arrow_submission_run_file(self):
        return op.join(self.log_dir, "submitted_arrow_jobs.txt")

    def _arrowed_bin_prefix(self, first, last):
        """Return $arrowed_dir/c{first}to{last}"""
        return self.arrowed_dir + "/c{first}to{last}".format(
            first=first, last=last)

    def fq_of_arrowed_bin(self, first, last):
        """Return $_arrowed_bin_prefix.arrowed.fastq
        this is arrowed fq output. Whenever this is changed, change
        IceArrowPostprocess2 accordingly.
        """
        return self._arrowed_bin_prefix(first, last) + ".arrowed.fastq"

    def list_of_expected_arrow_fq_files(self):
        """
        Used by IceArrowPostProcess2 to identify all the finished files.
        First looks at $root_dir/log/prepared_arrow_files.txt to identify all the c<i>to<j>.sh
        And figure out that the finished files should be c<i>to<j>.arrowed.fastq.

        Note: related methods include fq_of_arrowed_bin() and others,
              when changing one must remember to keep it consistent across.
        """
        def iter_script_to_get_fq(script_filename):
            for line in open(script_filename):
                # line might be like:
                # bash <arrow_dir>/c0to9.sh
                sh_file = line.strip().split()[-1]
                assert sh_file.endswith('.sh')
                yield sh_file[:-3] + '.arrowed.fastq'


        sge_ids = []
        submitted = {} # expected fq --> ("local" or SGE jobid, script used to get this)
        for line in open(self.arrow_submission_run_file):
            jobid, script = line.strip().split('\t')
            # read the script to see which c<i>to<j>.sh files are associated with this
            for fq in iter_script_to_get_fq(script):
                submitted[fq] = (jobid, script)
            if jobid!='local':
                sge_ids.append(jobid)

        return sge_ids, submitted


    @property
    def pattern_of_expected_fq_files(self):
        return "c*to*.arrowed.fastq"

    @property
    def unfinished_arrow_sh_files(self):
        """
        File that will contain the list of .sh files that were incomplete and need to re-run.
        """
        return op.join(self.log_dir, 'unfinished_arrow_files.txt')
