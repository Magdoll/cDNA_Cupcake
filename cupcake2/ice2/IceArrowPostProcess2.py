
"""
Post Arrow, pick up the high QV and low QV conesnsus isoforms.
"""

import os
import glob
import re
import logging
import os.path as op
from collections import defaultdict
from pickle import load
from time import sleep

from pbcore.io import FastqReader, FastqWriter, FastaWriter

from cupcake2.tofu2.ClusterOptions2 import IceArrowHQLQOptions2

from cupcake2.tofu2.ToFuOptions2 import  \
    add_cluster_summary_report_arguments, \
    add_ice_post_arrow_hq_lq_arguments2, \
    add_cluster_root_dir_as_positional_argument

from pbtranscript.Utils import phred_to_qv, as_contigset, \
    get_all_files_in_dir, ln, nfs_exists

from cupcake2.ice2.IceUtils2 import cid_with_annotation2
from cupcake2.ice2.IceFiles2 import IceFiles2
from cupcake2.ice2.__init__ import ICE_ARROW_PY


class IceArrowPostProcess2(IceFiles2):

    """check if arrow jobs are finished and arrow results are completed.
       If arrow jobs are completed, pick up high QV consensus isoforms.
       If arrow jobs are still running,
           * If quit_if_not_done is True, exit.
           * If quit_if_not_done is False, wait till arrow jobs are finished.

    """

    desc = "Post-polishing process, pick up high QV and low QV consensus isoforms."
    prog = "%s postprocess " % ICE_ARROW_PY

    def __init__(self, root_dir, ipq_opts,
                 quit_if_not_done=True,
                 summary_fn=None, report_fn=None,
                 no_log_f=False, make_dirs=True):
        super(IceArrowPostProcess2, self).__init__(
                prog_name="ice_arrow_postprocess2",
                root_dir=root_dir, no_log_f=no_log_f, make_dirs=make_dirs)
        self.quit_if_not_done = quit_if_not_done

        assert(type(ipq_opts) is IceArrowHQLQOptions2)
        self.ipq_opts = ipq_opts
        self.hq_isoforms_dataset = self.lq_isoforms_dataset = None
        self.hq_isoforms_fa = ipq_opts.hq_isoforms_fa
        if str(self.hq_isoforms_fa).endswith(".contigset.xml"):
            self.hq_isoforms_dataset = ipq_opts.hq_isoforms_fa
            self.hq_isoforms_fa = re.sub(".contigset.xml$", ".fasta",
                                         self.hq_isoforms_fa)
        self.hq_isoforms_fq = ipq_opts.hq_isoforms_fq
        self.lq_isoforms_fa = ipq_opts.lq_isoforms_fa
        if str(self.lq_isoforms_fa).endswith(".contigset.xml"):
            self.lq_isoforms_dataset = ipq_opts.lq_isoforms_fa
            self.lq_isoforms_fa = re.sub(".contigset.xml$", ".fasta",
                                         self.lq_isoforms_fa)
        self.lq_isoforms_fq = ipq_opts.lq_isoforms_fq
        # Arrow usually can't call accurate QV on both ends
        # because of very well + less coverage
        # Ignore QV of the first <100 bp> on 5' end
        self.qv_trim_5 = ipq_opts.qv_trim_5
        # Ignore QV of the last <30 bp> on the 3' end"""
        self.qv_trim_3 = ipq_opts.qv_trim_3
        # The max allowed average base error rate within
        # seq[qv_trim_5:-qv_trim_3]
        self.hq_arrow_min_accuracy = ipq_opts.hq_arrow_min_accuracy
        # The min number of supportive full-length reads for HQ
        self.hq_min_full_length_reads = ipq_opts.hq_min_full_length_reads

        self.fq_filenames = []

        self.report_fn = report_fn
        self.summary_fn = summary_fn

    def get_existing_binned_arrowed_fq(self):
        """
        Return all existing arrowed fq files for binned clusters.
        Each arrowed fq should be like c0to9.arrowed.fastq
        """
        return glob.glob(os.path.join(self.arrowed_dir, self.pattern_of_expected_fq_files))

    def validate_inputs(self):
        """Validate if logs and pickle for non-full-length reads exist."""
        errMsg = ""

        if not nfs_exists(self.nfl_all_pickle_fn):
            errMsg = "Pickle file {f} ".format(f=self.nfl_all_pickle_fn) + \
                     "which assigns all non-full-length reads to isoforms " + \
                     "does not exist. Please check 'ice_partial.py *' are " + \
                     "all done."
        elif not nfs_exists(self.final_pickle_fn):
            errMsg = "Pickle file {f} ".format(f=self.final_pickle_fn) + \
                     "which assigns full-length non-chimeric reads to " + \
                     "isoforms does not exist."
        elif not nfs_exists(self.arrow_submission_run_file):
            errMsg = "Log file {f}".format(f=self.arrow_submission_run_file) + \
                     " of all submitted arrow jobs does not exist."

        if errMsg != "":
            self.add_log(errMsg, level=logging.ERROR)
            raise IOError(errMsg)

    def check_arrow_jobs_completion(self):
        """Check whether arrow jobs are completed.
        submitted_arrow_jobs.txt should have format like:
        <job_id> \t ./arrowed/c0to10.sh

        Returns:
        "DONE" --- if all jobs are done and files are there return
        "FAILED" --- all jobs are done but some files incomplete ask if to resubmit
        "RUNNING" --- if not all jobs are done, just quit
        fq_filenames contains all the finished fastq files.
        """
        self.add_log("Checking if arrow jobs are completed.")
        bad_sh = []
        self.fq_filenames = []
        self.add_log("Submitted arrow jobs are at {f}:".
                     format(f=self.arrow_submission_run_file))

        # submitted = list of (SGE jobid or local, script file that is running)
        sge_jobids, submitted = self.list_of_expected_arrow_fq_files()

        done_flag = True
        running_jids = []
        # if one or more jobs were submitted through SGE,
        # go through qstat to see if anything is still running
        if len(sge_jobids):  # at least one job was run through SGE
            stuff = os.popen("qstat").read().strip().split('\n')
            assert stuff[0].startswith('job-ID')
            assert stuff[1].startswith('-------')
            # first two lines are header
            for x in stuff[2:]:
                job_id = x.split()[0]
                running_jids.append(job_id)
                if job_id in sge_jobids:
                    self.add_log("job {0} is still running.".format(job_id))
                    done_flag = False

        # now go through all the expected fastq files and check they exist
        for fq_filename,(job_id,sh_file) in submitted.items():
            if not nfs_exists(fq_filename) or \
                    os.stat(fq_filename).st_size == 0:
                if job_id in running_jids:  # still running, pass
                    done_flag = False
                else:
                    self.add_log("job {0} is completed but {1} is still empty!".
                                 format(job_id, fq_filename))
                    bad_sh.append(sh_file)
            else:
                self.add_log("job {0} is done".format(job_id))
                self.fq_filenames.append(fq_filename)

        if not done_flag:
            if len(bad_sh) == 0:
                return "RUNNING"
            else:
                # write the unfinished jobs to $unfinished_arrow_sh_files$
                f = open(self.unfinished_arrow_sh_files, 'w')
                f.write("\n".join(bad_sh) + '\n')
                f.close()
                self.add_log("Some jobs were incomplete! Please re-run all files listed in {1}.\n".format(\
                    len(bad_sh), f.name))
                return "FAILED"
        else:
            return "DONE"

    @property
    def arrowed_good_fa(self):
        """Return $root_dir/all_arrowed.hq.a_b_c.fasta"""
        return op.join(self.root_dir,
                       "all_arrowed_hq.{a}_{b}_{c}.fasta".format(
                           a=self.qv_trim_5,
                           b=self.qv_trim_3,
                           c=self.hq_arrow_min_accuracy))

    @property
    def arrowed_good_fq(self):
        """Return $root_dir/all_arrowed_hq.a_b_c.fastq"""
        return op.join(self.root_dir,
                       "all_arrowed_hq.{a}_{b}_{c}.fastq".format(
                           a=self.qv_trim_5,
                           b=self.qv_trim_3,
                           c=self.hq_arrow_min_accuracy))

    @property
    def arrowed_bad_fa(self):
        """Return $root_dir/all_arrowed_lq.fasta"""
        return op.join(self.root_dir, "all_arrowed_lq.fasta")

    @property
    def arrowed_bad_fq(self):
        """Return $root_dir/all_arrowed_lq.fastq"""
        return op.join(self.root_dir, "all_arrowed_lq.fastq")

    def pickup_best_clusters(self):
        """Pick up hiqh QV clusters."""
        self.add_log("Picking up the best clusters according to QVs from {fs}.".
                     format(fs=", ".join(self.fq_filenames)))
        a = load(open(self.final_pickle_fn))
        uc = a['uc']
        # check if the uc cids are integers
        uc_keys_are_int = type(list(uc.keys())[0]) is int

        polished = {} # cid --> FastqRecord

        for fq in self.fq_filenames:
            self.add_log("Looking at arrowed fq {f}".format(f=fq))
            for r in FastqReader(fq):
                # possible ID #1: c0|arrow (a single Ice2 directory)
                # possible ID #2: b112_c0|arrow (after collecting several Ice2 directory)
                cid = r.name.split('|')[0]
                if cid.endswith('_ref'):
                    cid = cid[:-4]
                i = cid.find('/')
                if i > 0:
                    cid = cid[:i]
                if uc_keys_are_int:
                    # only convert in the case where uc keys are integers (ex: is c10, but 10)
                    cid = int(cid[1:]) #becuz possible ID #2, dont convert to int
                polished[cid] = r


        expected_acc_dict = {} # cid --> expected accuracy (ex: 0.99)
        good = [] # contains all the cids that are HQ

        # calculate expected QV given 5'/3' trimming
        # for sequences that are shorter than the trimming, use the length itself
        for cid, r in polished.items():
            qv_len = max(len(r.quality), len(r.quality) - self.qv_trim_5 - self.qv_trim_3)
            q = [phred_to_qv(x) for x in r.quality]
            err_sum = sum(q[self.qv_trim_5: -self.qv_trim_3])
            expected_acc_dict[cid] = 1.0 - (err_sum / float(qv_len))
            if expected_acc_dict[cid] >= self.hq_arrow_min_accuracy and \
                len(uc[cid]) >= self.hq_min_full_length_reads :
                good.append(cid)

        partial_uc = load(open(self.nfl_all_pickle_fn))['partial_uc']
        partial_uc2 = defaultdict(lambda: [])
        partial_uc2.update(partial_uc)

        if self.report_fn is not None:
            self.write_report(report_fn=self.report_fn,
                              uc=uc, partial_uc=partial_uc2)

        self.add_log("Writing hiqh-quality isoforms to {f}|fq".
                     format(f=self.arrowed_good_fa))
        self.add_log("Writing low-quality isoforms to {f}|fq".
                     format(f=self.arrowed_bad_fa))
        with FastaWriter(self.arrowed_good_fa) as good_fa_writer, \
                FastaWriter(self.arrowed_bad_fa) as bad_fa_writer, \
                FastqWriter(self.arrowed_good_fq) as good_fq_writer, \
                FastqWriter(self.arrowed_bad_fq) as bad_fq_writer:
            for cid in polished:
                r = polished[cid]
                newname = "c{cid}/f{flnc_num}p{nfl_num}/{read_len}".\
                    format(cid=cid,
                           flnc_num=len(uc[cid]),
                           nfl_num=len(partial_uc2[cid]),
                           read_len=len(r.sequence))
                newname = cid_with_annotation2(newname, expected_acc=expected_acc_dict[cid])

                if cid in good:
                    self.add_log("processing arrowed cluster {c} --> good.".
                                 format(c=cid))
                    good_fa_writer.writeRecord(newname, r.sequence[:])
                    good_fq_writer.writeRecord(newname, r.sequence[:], r.quality)
                else:
                    self.add_log("processing arrowed cluster {c} --> bad.".
                                 format(c=cid))
                    bad_fa_writer.writeRecord(newname, r.sequence[:])
                    bad_fq_writer.writeRecord(newname, r.sequence[:], r.quality)

        self.add_log("-" * 60, level=logging.INFO)
        self.add_log("High-quality Arrowed consensus written " +
                     "to:\n{0}\n{1}".format(self.arrowed_good_fa,
                                            self.arrowed_good_fq),
                     level=logging.INFO)
        self.add_log("Low-quality Arrowed consensus written " +
                     "to:\n{0}\n{1}".format(self.arrowed_bad_fa,
                                            self.arrowed_bad_fq),
                     level=logging.INFO)
        self.add_log("-" * 60, level=logging.INFO)

    def cmd_str(self):
        """Return a cmd string ($ICE_ARROW_PY postprocess)."""
        return self._cmd_str(root_dir=self.root_dir,
                             ipq_opts=self.ipq_opts,
                             quit_if_not_done=self.quit_if_not_done,
                             summary_fn=self.summary_fn,
                             report_fn=self.report_fn)

    def _cmd_str(self, root_dir, ipq_opts, quit_if_not_done,
                 summary_fn, report_fn):
        """Return a cmd string ($ICE_ARROW_PY postprocess)."""
        cmd = self.prog + \
              "{d} ".format(d=root_dir) + \
              ipq_opts.cmd_str()
        if quit_if_not_done is True:
            cmd += "--quit_if_not_done "
        if summary_fn is not None:
            cmd += "--summary={f} ".format(f=summary_fn)
        if report_fn is not None:
            cmd += "--report={f} ".format(f=report_fn)
        return cmd

    def run(self):
        """
        Check all arrow jobs are running, failed or done. Write high-quality
        consensus and low-quality consensus to all_arrowed.hq|lq fasta|fastq
        """
        self.validate_inputs()

        job_stats = self.check_arrow_jobs_completion()
        self.add_log("Arrow job status: {s}".format(s=job_stats))

        if job_stats == 'DONE':
            pass # continue on below to process data
        elif job_stats == 'FAILED':
            self.add_log("Has incomplete jobs. Please re-run them.",
                         level=logging.ERROR)
            return -1
        elif job_stats == 'RUNNING':
            if self.quit_if_not_done:
                self.add_log("Jobs are still running. Please wait before running this script.")
                return 1
            else:
                while job_stats != "DONE":
                    self.add_log("Jobs are still running. Wait. Sleeping for 180 seconds.")
                    sleep(180)
                    job_stats = self.check_arrow_jobs_completion()
                    if job_stats == "DONE":
                        break
                    elif job_stats == "FAILED":
                        self.add_log("There are some failed jobs. Please check.",
                                     level=logging.ERROR)
                        return 1
                    elif job_stats == "RUNNING":
                        self.add_log("Jobs are still running. Wait. Sleeping for 180 seconds.",
                                     level=logging.INFO)
        else:
            msg = "Unable to recognize job_stats {s}".format(s=job_stats)
            self.add_log(msg, logging.ERROR)
            raise ValueError(msg)

        # at this point, all jobs must be done and all fastq files present.
        self.pickup_best_clusters()

        self.add_log("Creating polished high quality consensus isoforms.")
        if self.hq_isoforms_fa is not None:
            ln(self.arrowed_good_fa, self.hq_isoforms_fa)
        if self.hq_isoforms_fq is not None:
            ln(self.arrowed_good_fq, self.hq_isoforms_fq)

        self.add_log("Creating polished low quality consensus isoforms.")
        if self.lq_isoforms_fa is not None:
            ln(self.arrowed_bad_fa, self.lq_isoforms_fa)
        if self.lq_isoforms_fq is not None:
            ln(self.arrowed_bad_fq, self.lq_isoforms_fq)

        if self.hq_isoforms_dataset is not None:
            ds = as_contigset(self.hq_isoforms_fa, self.hq_isoforms_dataset)
        if self.lq_isoforms_dataset is not None:
            ds = as_contigset(self.lq_isoforms_fa, self.lq_isoforms_dataset)
        if self.summary_fn is not None:
            self.write_summary(summary_fn=self.summary_fn,
                               isoforms_fa=self.final_consensus_fa,
                               hq_fa=self.hq_isoforms_fa,
                               lq_fa=self.lq_isoforms_fa)

        self.close_log()


def add_ice_arrow_postprocess_arguments(arg_parser):
    """
    Add arguments for IceArrowPostProcess2 ($ICE_ARROW_PY postprocess).
    """
    arg_parser = add_cluster_root_dir_as_positional_argument(arg_parser)
    arg_parser = add_ice_post_arrow_hq_lq_arguments2(arg_parser)
    arg_parser = add_cluster_summary_report_arguments(arg_parser)

    arg_parser.add_argument("--quit_if_not_done",
                            default=False,
                            dest="quit_if_not_done",
                            action="store_true",
                            help="Quit if arrow jobs haven't been completed.")
    return arg_parser

