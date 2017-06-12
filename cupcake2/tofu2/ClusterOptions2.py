__author__ = 'etseng@pacb.com'

import os.path as op
import logging
import numpy as np
from pbtranscript.io.ContigSetReaderWrapper import ContigSetReaderWrapper

class SgeOptions2(object):

    """Define options to configure SGE."""

    def __init__(self, unique_id, use_sge=False, max_sge_jobs=40,
                 blasr_nproc=24, gcon_nproc=8, arrow_nproc=8,
                 sge_queue=None, sge_env_name="smp", qsub_extra=''):
        self.unique_id = unique_id
        self.use_sge = use_sge
        self.max_sge_jobs = max_sge_jobs
        self.blasr_nproc = blasr_nproc
        self.gcon_nproc = gcon_nproc
        self.arrow_nproc = arrow_nproc
        self.sge_queue = sge_queue
        self.sge_env_name = sge_env_name
        self.qsub_extra = qsub_extra  # whatever extra command customers want via qsub

    def __str__(self):
        return "unqiueID={i}\n".format(i=self.unique_id) + \
               "use_sge={u}\n".format(u=self.use_sge) + \
               "max_sge_jobs={j}\n".format(j=self.max_sge_jobs) + \
               "sge_queue={q}\n".format(q=self.sge_queue) + \
               "sge_env_name={s}\n".format(s=self.sge_env_name) + \
               "qsub_extra={e}\n".format(e=self.qsub_extra) + \
               "blasr_nproc={n}\n".format(n=self.blasr_nproc) + \
               "gcon_nproc={n}\n".format(n=self.gcon_nproc) + \
               "arrow_nproc={t}\n".format(t=self.arrow_nproc)

    def qsub_cmd(self, script, num_threads,
                 wait_before_exit=False, depend_on_jobs=None,
                 elog="/dev/null", olog="/dev/null",
                 is_script=True, jobid=None):
        """
        Return a cmd string to run script on sge.
        Parameters:
          script - script to run
          num_threads -  number of threads per sge node
          wait_before_exit - qsub wait for this job to complete before exiting
          depend_on_jobs - this job depends on a list of submitted jobs.
          elog - error log
          olog - output log
          is_script - if Yes, input is a script, else input is an executable.
        """
        # always execute the job from the current directory
        # -V specify all env vars
        # -S specify shell for the job
        # -pe {sge_env_name} <num_slots> specify parallel environment
        ret = "qsub -cwd -V -S /bin/bash "
        ret += "-pe {p} {n} ".format(p=self.sge_env_name,
                                     n=num_threads)
        if self.sge_queue is not None:
            ret += "-q {q} ".format(q=self.sge_queue)

        if self.qsub_extra != '':
            ret += "{e} ".format(e=self.qsub_extra)

        if jobid is not None:
            ret += "-N {j} ".format(j=jobid)

        if wait_before_exit is True:
            ret += "-sync y "

        if depend_on_jobs is not None:
            ret += "-hold_jid {jids} ".format(jids=",".join(depend_on_jobs))

        ret += "-e {e} -o {o} ".format(e=elog, o=olog)

        if not is_script:
            ret += "-b y "

        ret += "{script}".format(script=script)
        return ret

    def cmd_str(self, show_blasr_nproc=False, show_gcon_nproc=False,
                show_arrow_nproc=False, show_sge_queue=False,
                show_sge_env_name=False, show_qsub_extra=False):
        """Return a cmd string."""
        cmd = ""
        if self.use_sge is True:
            cmd += "--use_sge "
            cmd += "--max_sge_jobs={n} ".format(n=self.max_sge_jobs)
            cmd += "--unique_id={n} ".format(n=self.unique_id)
        if show_blasr_nproc is True and self.blasr_nproc is not None:
            cmd += "--blasr_nproc={n} ".format(n=self.blasr_nproc)
        if show_gcon_nproc is True and self.gcon_nproc is not None:
            cmd += "--gcon_nproc={n} ".format(n=self.gcon_nproc)
        if show_arrow_nproc is True and self.arrow_nproc is not None:
            cmd += "--arrow_nproc={n} ".format(n=self.arrow_nproc)
        if show_sge_env_name is True:
            cmd += "--sge_env_name={0}".format(self.sge_env_name)
        if show_sge_queue is True and self.sge_queue is not None:
            cmd += "--sge_queue={0}".format(self.sge_queue)
        if show_qsub_extra and self.qsub_extra!='':
            cmd += "--qsub_extra={0}".format(self.qsub_extra)
        return cmd


class IceOptions2(object):
    """
    Define ICE related options.

    aligner_choice: daligner, blasr
    """

    def __init__(self, cDNA_size="under1k", flnc_reads_per_split=20000,
                 ece_penalty=1, ece_min_len=20, bestn=100,
                 max_missed_start=400, max_missed_end=100,
                 min_match_len=50,
                 quiver=False,
                 use_finer_qv=False, targeted_isoseq=False,
                 nfl_reads_per_split=30000,
                 num_clusters_per_bin=100,
                 aligner_choice='daligner'):

        self.cDNA_size = str(cDNA_size)

        self.low_cDNA_size = None
        self.high_cDNA_size = None

        self.aligner_choice = aligner_choice

        # (user-set) maximum allowed missed/start to be considered an "isoform"
        # setting this allows for 5' degradation and some slight 3' differences
        # recommended: 400 bp for 5', 50 bp for 3'
        self.max_missed_start = int(max_missed_start)
        self.max_missed_end = int(max_missed_end)

        # maximum allowed missed start/end to be considered "fully" aligned by an aligner
        # must be smaller than max_missed_start/end, recommend to keep it tight
        # setting it too high (like 200 bp) may mean accidentally clustering diff isoforms
        # that have different exonic structure in the first 200 bp
        self.full_missed_start = 50
        self.full_missed_end = 30

        self.ece_penalty = int(ece_penalty)
        self.ece_min_len = int(ece_min_len)
        self.bestn = int(bestn)
        self.quiver = quiver
        self.targeted_isoseq = targeted_isoseq
        # flnc reads per split
        self.flnc_reads_per_split = int(flnc_reads_per_split)
        # nfl reads per split
        self.nfl_reads_per_split = int(nfl_reads_per_split)
        # whether or not to use finer qv from ccs.h5 (including
        # deletionQV, insertionQV, mergeQV) or a single qv from fastq
        self.use_finer_qv = use_finer_qv
        # Put every 100 clusters in to a bin for quiver
        self.num_clusters_per_bin = num_clusters_per_bin

        self.min_match_len = min_match_len # by default, min_match_len (used by aligners) is 50 bp

    @classmethod
    def cDNA_sizeBins(cls):
        """Return cDNA size bins."""
        return ("under1k", "between1k2k", "between2k3k", "3to5k", "above5k")

    def detect_cDNA_size(self, fasta_filename):
        """
        Auto detection cDNA_size adjust maxScore and sensitive mode accordingly.
        Parameters:
          fasta_filename - use 90% percentile read length as cDNA size.
        """
        if not (fasta_filename.lower().endswith(".fasta") or
                fasta_filename.lower().endswith(".fa") or
                fasta_filename.lower().endswith(".xml")):
            raise ValueError("%s must take a FASTA ContigSet file as input."
                             % fasta_filename)

        if op.exists(self._config_filename(fasta_filename)):
            self._read_config(fasta_filename=fasta_filename)
        else:
            self._write_config(fasta_filename=fasta_filename)

        if self.low_cDNA_size <= 1000:
            self.cDNA_size = "under1k"
        elif  self.low_cDNA_size <= 2000:
            self.cDNA_size = "between1k2k"
        elif self.low_cDNA_size <= 3000:
            self.cDNA_size = "between2k3k"
        elif self.low_cDNA_size <= 5000:
            self.cDNA_size = "3to5k"
        else:
            self.cDNA_size = "above5k"

        logging.info("Auto-detection of cDNA size set to {n} (maxScore: {s})"
                     .format(n=self.cDNA_size, s=self.maxScore))

    @property
    def sensitive_mode(self):
        """Return whether or not use sensitive mode.

        DALIGNER uses sensitive mode when most reads are >= 6K
        """
        return self.low_cDNA_size is None or \
               self.low_cDNA_size >= 6000

    def _config_filename(self, fasta_filename):
        """ return config filename = fasta_filename.sensitive.config """
        return "%s.sensitive.config" % fasta_filename

    def _write_config(self, fasta_filename):
        """Write daligner sensitive config to fasta_filename.sensitive.config."""
        lens = [len(r.sequence) for r in ContigSetReaderWrapper(fasta_filename)]
        self.low_cDNA_size, self.high_cDNA_size = 0, 0
        if len(lens) == 1:
            self.low_cDNA_size, self.high_cDNA_size = lens[0], lens[0]
        if len(lens) >= 2:
            self.low_cDNA_size  = int(np.percentile(lens, 10))
            self.high_cDNA_size = int(np.percentile(lens, 90))

        try:
            with open(fasta_filename+'.sensitive.config', 'w') as f:
                f.write("sensitive={s}\n".format(s=self.sensitive_mode))
                f.write("low={l}\n".format(l=self.low_cDNA_size))
                f.write("high={h}\n".format(h=self.high_cDNA_size))
        except IOError:
            pass # it's OK not to have write permission

    def _read_config(self, fasta_filename):
        """Read from sensitive.config file."""
        if not op.exists(self._config_filename(fasta_filename)):
            raise IOError("Could not find config file: %s." % fasta_filename)
        else:
            with open(self._config_filename(fasta_filename)) as f:
                a, b = f.readline().strip().split('=')
                assert a == 'sensitive'
                #flag = (b == 'True')
                a, b = f.readline().strip().split('=')
                assert a == 'low'
                self.low_cDNA_size = int(b)
                a, b = f.readline().strip().split('=')
                assert a == 'high'
                self.high_cDNA_size = int(b)

    @property
    def maxScore(self):
        """Return maximum blasr score according to estimated cDNA size."""
        if self.cDNA_size not in IceOptions2.cDNA_sizeBins():
            raise ValueError("Invalid cDNA size: {cs}".
                             format(cs=self.cDNA_size))
        d = {"under1k": -1000, "between1k2k": -2000, "between2k3k": -3000,
             "3to5k": -5000, "above5k": -8000}
        return d[self.cDNA_size]

    def cmd_str(self):
        return self.__str__()

    def __str__(self):
        return "cDNA_size={sz}\n".format(sz=self.cDNA_size) + \
               "" if self.low_cDNA_size is None else \
               "low={l}\nsensitive={s}\n".format(l=self.low_cDNA_size,
                                                 s=self.sensitive_mode) + \
               "" if self.high_cDNA_size is None else \
               "high={h}\n".format(h=self.high_cDNA_size) + \
               "maxScore={ms}\n".format(ms=self.maxScore) + \
               "ece_penalty={ep}\n".format(ep=self.ece_penalty) + \
               "ece_min_len={eml}\n".format(eml=self.ece_min_len) + \
               "bestn={bsn}\n".format(bsn=self.bestn) + \
               "aligner_choice={0}\n".format(self.aligner_choice) + \
               "max_missed_start={0}\n".format(self.max_missed_start) + \
               "max_missed_end={0}\n".format(self.max_missed_end) + \
               "quiver={quiver}\n".format(quiver=self.quiver) + \
               "targeted_isoseq={t}\n".format(t=self.targeted_isoseq) + \
               "flnc_reads_per_split={n}\n".format(n=self.flnc_reads_per_split) + \
               "use_finer_qv={qv}\n".format(qv=self.use_finer_qv) + \
               "nfl_reads_per_split={n}\n".format(n=self.nfl_reads_per_split) + \
               "num_clusters_per_bin={n}\n".format(n=self.num_clusters_per_bin)


class IceArrowHQLQOptions2(object):

    """
    Define HQ/LQ isoforms related options

    Liz Note: we are changing `hq_min_full_length_reads` potentially back to min_FL>=1 again!
              evaluating how we can increase observed isoform # while still removing artifacts
    """

    def __init__(self, qv_trim_5=100, qv_trim_3=30, hq_arrow_min_accuracy=0.99,
                 hq_isoforms_fa=None, hq_isoforms_fq=None,
                 lq_isoforms_fa=None, lq_isoforms_fq=None,
                 hq_min_full_length_reads=1):
        # Ignore QV of n bases in the 5' end
        self.qv_trim_5 = int(qv_trim_5)
        # Ignore QV of n bases in the 3' end
        self.qv_trim_3 = int(qv_trim_3)
        # Minimum allowed arrow accuracy to mark an isoform as HQ
        self.hq_arrow_min_accuracy = float(hq_arrow_min_accuracy)

        assert 0 < self.hq_arrow_min_accuracy <= 1

        self.hq_isoforms_fa = hq_isoforms_fa
        self.hq_isoforms_fq = hq_isoforms_fq
        self.lq_isoforms_fa = lq_isoforms_fa
        self.lq_isoforms_fq = lq_isoforms_fq
        self.hq_min_full_length_reads = hq_min_full_length_reads

    def __str__(self):
        return "qv_trim_5={n}\n".format(n=self.qv_trim_5) + \
               "qv_trim_3={n}\n".format(n=self.qv_trim_3) + \
               "hq_arrow_min_accuracy={n}\n".\
               format(n=self.hq_arrow_min_accuracy) + \
               "hq_min_full_length_reads={n}\n".\
               format(n=self.hq_min_full_length_reads) + \
               "HQ isoforms fasta={fa}\n".format(fa=self.hq_isoforms_fa) + \
               "HQ isoforms fastq={fq}\n".format(fq=self.hq_isoforms_fq) + \
               "LQ isoforms fasta={fa}\n".format(fa=self.lq_isoforms_fa) + \
               "LQ isoforms fastq={fq}\n".format(fq=self.lq_isoforms_fq)

    def cmd_str(self):
        """Return a cmd string."""
        cmd = "--hq_arrow_min_accuracy={n} ".\
                format(n=self.hq_arrow_min_accuracy) + \
              "--qv_trim_5={n} ".format(n=self.qv_trim_5) + \
              "--qv_trim_3={n} ".format(n=self.qv_trim_3)
        if self.hq_isoforms_fa is not None:
            cmd += "--hq_isoforms_fa={f} ".format(f=self.hq_isoforms_fa)
        if self.hq_isoforms_fq is not None:
            cmd += "--hq_isoforms_fq={f} ".format(f=self.hq_isoforms_fq)
        if self.lq_isoforms_fa is not None:
            cmd += "--lq_isoforms_fa={f} ".format(f=self.lq_isoforms_fa)
        if self.lq_isoforms_fq is not None:
            cmd += "--lq_isoforms_fq={f} ".format(f=self.lq_isoforms_fq)
        return cmd