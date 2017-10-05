#!/usr/bin/env python

"""
Polish consensus isoforms created by ICE, using Arrow for Sequel data.

Steps:
1. collect all the <cid>s that need to be polished by arrow (ex: c0, c2, c10....c200)
2. determine which <cid>s go into each job file (ex: c0to100, c301to400), write out as text
3. make arrowed/c<i>to<j>.arrowed.sh for each job file, write out as text
4a. if local, run through each job sequentially
4b. if using SGE, partition the total jobs from (3) by max_sge_jobs
"""

import os.path as op
import logging
import shutil
import cPickle
import json
from math import ceil
from collections import defaultdict

from pbtranscript.ClusterOptions import IceQuiverOptions

from pbtranscript.Utils import mkdir, real_upath, nfs_exists, \
    get_files_from_file_or_fofn, guess_file_format, FILE_FORMATS, \
    use_samtools_v_1_3_1
from pbtranscript.ice.IceUtils import get_the_only_fasta_record, \
    is_blank_sam, concat_sam, blasr_for_quiver, trim_subreads_and_write, \
    is_blank_bam, concat_bam

from cupcake2.tofu2.ToFuOptions2 import add_fofn_arguments, \
    add_sge_arguments, add_cluster_root_dir_as_positional_argument, BaseConstants
from cupcake2.ice2.IceFiles2 import IceFiles2


from pbtranscript.io import MetaSubreadFastaReader, BamCollection, \
    FastaRandomReader
from pbcore.io import FastaWriter


class IceArrow2(IceFiles2):
    """
    IceArrow2
    """

    desc = "After assigning all non-full-length reads to unpolished " + \
           "consensus isoforms created by ICE, polish these consensus " + \
           "isoforms, using arrow for Sequel data."

    def __init__(self, root_dir, subread_xml, sge_opts,
                 tmp_dir=None, prog_name=None):
        # Initialize super class IceFiles.
        prog_name = "IceArrow2" if prog_name is None else prog_name
        super(IceArrow2, self).__init__(prog_name=prog_name,
                                        root_dir=root_dir, subread_xml=subread_xml,
                                        tmp_dir=tmp_dir)
        self.sge_opts = sge_opts
        self.use_samtools_v_1_3_1 = use_samtools_v_1_3_1()

    def validate_inputs(self):
        """Validate input fofns, and root_dir, log_dir, tmp_dir,
        create arrowed_dir and arrowed_log_dir"""
        self.add_log("Validating inputs.")

        # Create directories: root_dir/quivered and root_dir/log_dir/quivered
        try:
            mkdir(self.arrowed_dir)
            mkdir(self.arrowed_log_dir)
        except OSError:
            # Multiple ice_arrow_i jobs may run at the same time and try to
            # mkdir, race condition may happen, so ignore OSError here.
            pass

        errMsg = ""

        if not nfs_exists(self.log_dir) or not op.isdir(self.log_dir):
            errMsg = "Log dir {l} is not an existing directory.".\
                format(l=self.log_dir)
        elif self.subread_xml is None:
            errMsg = "Please specify subreads XML (e.g., --subread_xml=<movie>.subreadset.xml)."
        elif not nfs_exists(self.subread_xml):
            errMsg = "Specified subreads file (subread_xml={f}) does not exist.".format(f=self.subread_xml)
        elif guess_file_format(self.subread_xml) is not FILE_FORMATS.BAM:
            errMsg = "Invalid subreads XML file: {0}!".format(self.subread_xml)
        elif not nfs_exists(self.nfl_all_pickle_fn):
            #"output/map_noFL/noFL.ALL.partial_uc.pickle"):
            errMsg = "Pickle file {f} ".format(f=self.nfl_all_pickle_fn) + \
                     "which assigns all non-full-length reads to isoforms " + \
                     "does not exist. Please check 'run_IcePartials2.py *' are " + \
                     "all done."
        elif not nfs_exists(self.final_pickle_fn):
            errMsg = "Pickle file {f} ".format(f=self.final_pickle_fn) + \
                     "which assigns full-length non-chimeric reads to " + \
                     "isoforms does not exist."

        if errMsg != "":
            self.add_log(errMsg, level=logging.ERROR)
            raise IOError(errMsg)


    def sam_of_arrowed_bin(self, first, last):
        """Return $_arrowed_bin_prefix.sam"""
        return self._arrowed_bin_prefix(first, last) + ".sam"

    def bam_of_arrowed_bin(self, first, last, is_sorted=False):
        """
        Return $_arrowed_bin_prefix.unsorted.bam if not sorted;
        return $_arrowed_bin_prefix.bam if sorted.
        """
        if not is_sorted:
            return self._arrowed_bin_prefix(first, last) + ".unsorted.bam"
        else:
            return self._arrowed_bin_prefix(first, last) + ".bam"

    def ref_fa_of_arrowed_bin(self, first, last):
        """Return $_arrowed_bin_prefix.ref.fasta
        this is reference fasta for quiver to use as input.
        """
        return self._arrowed_bin_prefix(first, last) + ".ref.fasta"



    def script_of_arrowed_bin(self, first, last):
        """Return $_arrowed_bin_prefix.sh"""
        return self._arrowed_bin_prefix(first, last) + ".sh"

    def prepare_clusters_per_arrow_file(self, cids):
        """
        Return a list of [(first_index, last_index)] that goes in each chunk
        ex: [(1, 100), (101, 200), (201, 300) etc...]
        meaning we will have c1to100.arrowed.sh, c156to301.arrowed.sh, ...
        """
        n = BaseConstants.HQ_ARROW_CIDS_PER_FILE
        m = len(cids)
        return [(i,min(m,i+n)) for i in xrange(0, m, n)]

    def reconstruct_ref_fa_for_clusters_in_bin(self, cids, refs):
        """
        Reconstruct ref_fa of the cluster in the new tmp_dir
        e.g.,
            self.g_consensus_ref_fa_of_cluster(cid)

        Liz: new cids after ice2 collection is b<bin>_c<cid>
        refs --- dict{int(cid): ref_fa of cluster(cid)}
        """
        # Check existence when first time it is read.
        if not nfs_exists(self.final_consensus_fa):
            raise IOError("Final consensus FASTA file {f}".format(
                f=self.final_consensus_fa) + "does not exist.")

        print "Reconstructing g consensus files for clusters {0}, {1} in {2}".format(cids[0], cids[-1], self.tmp_dir)
        self.add_log("Reconstructing g consensus files for clusters {0}, {1} in {2}".format(cids[0], cids[-1], self.tmp_dir))

        final_consensus_d = FastaRandomReader(self.final_consensus_fa)
        for ref_id in final_consensus_d.d.keys():
            # Liz: this is no longer valid for the Ice2 cids #cid = int(ref_id.split('/')[0].replace('c', ''))
            cid = ref_id
            if cid in cids:
                _dir = self.cluster_dir_for_reconstructed_ref(cid)
                mkdir(_dir)
                ref_fa = op.join(_dir, op.basename(refs[cid]))
                refs[cid] = ref_fa
                with FastaWriter(ref_fa) as writer:
                    self.add_log("Writing ref_fa %s" % refs[cid])
                    writer.writeRecord(ref_id,
                                       final_consensus_d[ref_id].sequence[:])

        self.add_log("Reconstruct of g consensus files completed.",
                     level=logging.INFO)

    def create_raw_files_for_clusters_in_bin(self, cids, d, uc, partial_uc, refs):
        """
        Create raw subreads bam files for clusters in cids.
        For each cluster k in cids,
        * Collect raw subreads of zmws associated with cluster k
          in either uc or partial_uc.

        cids --- cluster ids
        d --- BamCollection
        uc --- uc[k] returns fl ccs reads associated with cluster k
        partial_uc --- partial_uc[k] returns nfl ccs reads associated with cluster k
        """
        file_func = self.raw_bam_of_cluster2

        for k in cids:  # for each cluster k
            # write cluster k's associated raw subreads to raw_fa
            # Trim both ends of subreads (which contain primers and polyAs)
            trim_subreads_and_write(reader=d,
                                    in_seqids=uc[k] + partial_uc[k],
                                    out_file=file_func(refs[k]),
                                    trim_len=IceQuiverOptions.trim_subread_flank_len,
                                    min_len=IceQuiverOptions.min_trimmed_subread_len,
                                    ignore_keyerror=True,
                                    bam=True)

    def create_sams_for_clusters_in_bin(self, cids, refs):
        """
        Create sam files for clusters in cids.
        For each cluster k in cids,
        * Call blasr to align its associated subreads to its consensus
          sequence as reference.

        cids --- cluster ids
        refs --- refs[k] -> consensus seq of cluster k

        This function has to be called after raw_fa_of_cluster for clusters
        in cids are created.

        """
        raw_file_func  = self.raw_bam_of_cluster2
        out_file_func = self.bam_of_cluster2

        for k in cids:  # for each cluster k
            # $root_dir/tmp/?/c{k}/in.raw_with_partial.fasta
            raw_fn = raw_file_func(refs[k])
            out_fn = out_file_func(refs[k])

            if not op.exists(raw_fn):
                raise IOError("{f} does not exist. ".format(f=raw_fn) +
                              "Please check raw subreads of this bin is created.")
            blasr_for_quiver(
                query_fn=raw_fn,
                ref_fasta=refs[k],
                out_fn=out_fn,
                bam=True,
                run_cmd=True,
                blasr_nproc=self.sge_opts.blasr_nproc)


    def concat_valid_sams_and_refs_for_bin(self, cids, refs):
        """
        Concatenate sam files and reference sequences of all valid clusters
        in bin to create a big sam and a big ref.
        A cluster is not valid if (1) or (2)
            (1) identical sequences already exists in another cluster in the same file
                (rare, but happens)
            (2) the alignment is empty (also rare, but happens)
        Return valid_cids, a list of valid cluster ids
        """
        first, last = cids[0], cids[-1]
        bin_ref_fa = self.ref_fa_of_arrowed_bin(first, last)
        bin_sam_file = self.bam_of_arrowed_bin(first, last)
        file_func  = self.bam_of_cluster2
        is_blank_file = is_blank_bam
        concat_sambam = concat_bam

        self.add_log("Concatenating reference files between " +
                     "{first} and {last}.".format(first=first, last=last))
        valid_sam_files = []
        valid_cids = []
        seqs_seen = {}
        with open(bin_ref_fa, 'w') as bin_ref_fa_writer:
            for cid in cids:
                fname = file_func(refs[cid])
                if not is_blank_file(fname):
                    ref_rec = get_the_only_fasta_record(refs[cid])
                    name = ref_rec.name.strip()
                    seq = ref_rec.sequence.strip()
                    if seq not in seqs_seen:
                        valid_sam_files.append(fname)
                        valid_cids.append(cid)
                        seqs_seen[seq] = cid
                        # concate valid ref files, avoid 'cat ...' hundreds
                        # or even thousands of files due to linux cmd line
                        # length limits
                        bin_ref_fa_writer.write(">{0}\n{1}\n".
                                                format(name, seq))
                    else:
                        self.add_log("ignoring {0} because identical " +
                                     "sequence!".format(cid))
                else:
                    self.add_log(
                        "ignoring {0} because no alignments!".format(cid))

        if len(valid_sam_files) == 0:
            self.add_log("No alignments were found for clusters between " +
                         "{first} and {last}.".format(first=first, last=last),
                         level=logging.WARNING)
            assert(len(valid_cids) == 0)
        else:
            self.add_log("Concatenating sam files between " +
                         "{first} and {last}.".format(first=first, last=last))
            # concat valid sam files
            concat_sambam(valid_sam_files, bin_sam_file)
            self.add_log("Concatenation done")

        return valid_cids

    def arrow_cmds_for_bin(self, cids):
        """
        Return a list of quiver related cmds. Input format must be BAM.
        """
        first, last = cids[0], cids[-1]
        self.add_log("Creating arrow cmds for c{first} to c{last}".
                     format(first=first, last=last))

        bin_ref_fa = self.ref_fa_of_arrowed_bin(first, last)
        bin_fq = self.fq_of_arrowed_bin(first, last)

        bin_unsorted_bam_file = self.bam_of_arrowed_bin(first, last, is_sorted=False)
        bin_bam_file = self.bam_of_arrowed_bin(first, last, is_sorted=True)
        bin_bam_prefix = self._arrowed_bin_prefix(first, last)

        cmds = []
        if not self.use_samtools_v_1_3_1:
            # SA2.*, SA3.0, SA3.1 and SA3.2 use v0.1.19
            cmds.append("samtools sort {f} {d}".format(
                f=real_upath(bin_unsorted_bam_file),
                d=real_upath(bin_bam_prefix)))
        else:
            # SA3.3 and up use v1.3.1
            cmds.append("samtools sort {f} -o {d}.bam".format(
                f=real_upath(bin_unsorted_bam_file),
                d=real_upath(bin_bam_prefix)))

        cmds.append("samtools index {f}".format(f=real_upath(bin_bam_file)))
        cmds.append("samtools faidx {ref}".format(ref=real_upath(bin_ref_fa)))
        cmds.append("pbindex {f}".format(f=real_upath(bin_bam_file)))
        cmds.append("variantCaller --maskRadius 3 -x 1 --minAccuracy 0 --algorithm=best " +
                    "{f} ".format(f=real_upath(bin_bam_file)) +
                    "--verbose -j{n} ".format(n=self.sge_opts.arrow_nproc) +
                    "--referenceFilename={ref} ".format(ref=real_upath(bin_ref_fa)) +
                    "-o {fq}".format(fq=real_upath(bin_fq)))
        return cmds

    def create_arrow_sh_for_bin(self, cids, cmds):
        """
        Write cmds to a bash script, e.g., arrowed/c{}to{}.sh,
        return script file path.
        """
        first, last = cids[0], cids[-1]
        bin_sh = self.script_of_arrowed_bin(first, last)
        self.add_log("Creating arrow bash script {f} for c{first} to c{last}.".
                     format(f=bin_sh, first=first, last=last))
        with open(bin_sh, 'w') as f:
            f.write("#!/bin/bash\n")
            f.write("\n".join(cmds))
        return bin_sh

    def create_submission_jobs(self, arrow_sh_scripts):
        """
        arrow_sh_scripts: the full list of arrowed/c{}to{}.sh we want to run.

        if using SGE, break into <max_sge_jobs>.
        if not using SGE, write into a big file where each is just:

        bash c{}to{}.sh
        bash c{}to{}.sh
        ...
        """
        files = []
        total_jobs = min(len(arrow_sh_scripts), self.sge_opts.max_sge_jobs if self.sge_opts.use_sge else 1)
        script_per_job = len(arrow_sh_scripts) / total_jobs + (1 if len(arrow_sh_scripts)%total_jobs > 0 else 0)

        for i in xrange(total_jobs):
            with open(self.arrow_submission_file(i, total_jobs), 'w') as f:
                for j in xrange(i*script_per_job, min(len(arrow_sh_scripts), (i+1)*script_per_job)):
                    f.write("bash {0}\n".format(arrow_sh_scripts[j]))
                files.append(f.name)

        return files


    def submit_jobs_local_or_remote(self, files_to_run):
        """
        Run jobs either locally or through SGE.
        Return a list of [(sge_job_id, filename)], which
        is also written to log/submitted_arrow_jobs.txt
        """
        flag_run_locally = (self.sge_opts.use_sge is not True) or (self.sge_opts.max_sge_jobs==0)
        if flag_run_locally:
            self.add_log("Files to submit locally: {0}\n".format(",".join(files_to_run)))
        else:
            self.add_log("Files to submit through SGE: {0}\n".format(",".join(files_to_run)))


        submit_f = open(self.arrow_submission_run_file, 'w')

        submitted = []
        for file in files_to_run:
            elog = op.join(self.arrowed_log_dir, op.basename(file) + ".elog")
            olog = op.join(self.arrowed_log_dir, op.basename(file) + ".olog")

            if flag_run_locally:
                cmd = "bash {f}".format(f=real_upath(file))
                self.run_cmd_and_log(cmd, olog=olog, elog=elog,
                                     description="Failed to run Arrow")
                submitted.append(("local", file))
                submit_f.write("{0}\t{1}\n".format("local", file))
            else:
                jid = "ice_arrow_{unique_id}_{name}".format(
                    unique_id=self.sge_opts.unique_id,
                    name=op.basename(file))
                qsub_cmd = self.sge_opts.qsub_cmd(script=file,
                                                  num_threads=self.sge_opts.arrow_nproc,
                                                  wait_before_exit=False,
                                                  depend_on_jobs=None,
                                                  elog=elog,
                                                  olog=olog,
                                                  is_script=True,
                                                  jobid=jid)
                job_id = self.qsub_cmd_and_log(qsub_cmd)
                submitted.append((job_id, file))
                submit_f.write("{0}\t{1}\n".format(job_id, file))
        submit_f.close()
        return submitted

    def create_a_arrow_bin(self, cids, d, uc, partial_uc, refs):
        """Put clusters in cids together into a bin. In order to polish
        consensus of clusters in the bin, prepare inputs and create a quiver
        bash script to run later.

        (1) For each cluster k in cids, obtain subreads of all zmws
            belonging to this cluster, and save in raw_fa_of_cluster(k)
        (2) For each cluster k in cids, call blasr to align raw_fa_of_cluster to
            its consensus sequence and create sam_of_cluster(k).
        (3) Concat all sam files of `valid` clusters to sam_of_quivered_bin, and
            concat ref seqs of all `valid` clusters to ref_fa_of_arrowed_bin
        (4) Make commands including
                samtoh5, loadPulses, cmph5tools.py, loadChemistry, ..., quiver
            in order to convert sam_of_quivered_bin to cmph5_of_quivered_bin.
            Write these commands to script_of_arrowed_bin
              * qsub all jobs later when scripts of all quivered bins are done.
              * or execute scripts sequentially on local machine
        """
        if not isinstance(d, BamCollection):
            raise TypeError("%s.create_a_arrow_bin, does not support %s" %
                            (self.__class__.__name__, type(d)))

        self.add_log("Creating a arrow job bin for clusters "
                     "[%s, %s]" % (cids[0], cids[-1]), level=logging.INFO)

        # For each cluster in bin, create its raw subreads fasta file.
        self.create_raw_files_for_clusters_in_bin(cids=cids, d=d, uc=uc,
                                                  partial_uc=partial_uc, refs=refs)

        # For each cluster in bin, align its raw subreads to ref to build a sam
        self.create_sams_for_clusters_in_bin(cids=cids, refs=refs)

        # Concatenate sam | ref files of 'valid' clusters in this bin to create
        # a big sam | ref file.
        valid_cids = self.concat_valid_sams_and_refs_for_bin(cids=cids,
                                                             refs=refs)

        # quiver cmds for this bin
        if len(valid_cids) != 0:
            cmds = self.arrow_cmds_for_bin(cids=cids)
        else:
            cmds = ["echo no valid clusters in this bin, skip..."]

        # Write quiver cmds for this bin to $root_dir/quivered/c{}_{}.sh
        return self.create_arrow_sh_for_bin(cids=cids, cmds=cmds)


    def create_arrows_bins_no_submit(self, d, uc, partial_uc, refs, cids_todo):
        """
        Create arrow bins for cids in <cids_todo>. Handle missing references, etc.
        Create/Write the jobs but DO NOT submit.
        """
        # Liz: I'm commenting this out because the "refs" from the pickle should be accurate
        # plus the new cids after ice2 collection is b<bin>_c<cid>
        # Update refs
        #new_refs = {cid: op.join(self.cluster_dir(cid), op.basename(refs[cid])) for cid in cids_todo}
        #refs = new_refs

        #print "create_arrows_bins_no_submit calld for {0}-{1}, {2} files".format(cids_todo[0],cids_todo[-1], len(cids_todo))
        # Reconstruct refs if not exist.
        cids_missing_refs = filter(lambda x: not op.exists(refs[x]), cids_todo)
        #print "{0} missing refs".format(len(cids_missing_refs))
        if len(cids_missing_refs) > 0:
            self.reconstruct_ref_fa_for_clusters_in_bin(cids=cids_missing_refs, refs=refs)

        return self.create_a_arrow_bin(cids_todo, d, uc, partial_uc, refs)


    @property
    def report_fn(self):
        """Return a csv report with cluster_id, read_id, read_type."""
        return op.join(self.out_dir, "cluster_report.FL_nonFL.csv")

    def load_pickles(self):
        """Load uc and refs from final_pickle_fn, load partial uc from
        nfl_all_pickle_fn, return (uc, partial_uc. refs).
        """
        def _load_pickle(fn):
            """Load *.json or *.pickle file."""
            with open(fn) as f:
                if fn.endswith(".json"):
                    return json.loads(f.read())
                else:
                    return cPickle.load(f)
        self.add_log("Loading uc from {f}.".format(f=self.final_pickle_fn))
        a = _load_pickle(self.final_pickle_fn)
        uc = a['uc']
        refs = a['refs']

        self.add_log("Loading partial uc from {f}.".
                     format(f=self.nfl_all_pickle_fn))
        partial_uc = _load_pickle(self.nfl_all_pickle_fn)['partial_uc']
        partial_uc2 = defaultdict(lambda: [])
        partial_uc2.update(partial_uc)
        return (uc, partial_uc2, refs)

    def index_input_subreads(self):
        """Index input subreads in self.subread_xml
        """
        msg = "Indexing files in %s, please wait." % self.subread_xml
        self.add_log(msg)
        d = BamCollection(self.subread_xml)
        self.add_log("File indexing done.")
        return d


    def run(self):
        """
        Run Arrow to polish all consensus isoforms predicted by ICE.
        1. load FL + nFL association information from pickle
        2. index subreads XML file
        """
        # Validate inputs
        self.validate_inputs()

        # read pickle and write out cluster_report.csv
        uc, partial_uc, refs = self.load_pickles()
        self.write_report(report_fn=self.report_fn,
                              uc=uc, partial_uc=partial_uc)

        # good = [x for x in uc if len(uc[x]) > 1 or len(partial_uc2[x]) >= 10]
        # bug 24984, call quiver on everything, no selection is needed.
        cid_keys = sorted([x for x in uc])  # sort cluster ids
        cid_index_chunks = self.prepare_clusters_per_arrow_file(cid_keys)

        # write out to a text file so we know the cid chunks if jobs fail at this point
        with open(self.prepare_arrow_files, 'w') as f:
            for first_index,last_index in cid_index_chunks:
                f.write(self.script_of_arrowed_bin(cid_keys[first_index], cid_keys[last_index-1]) + '\n')

        # Index input <subread_xml> and make the actual arrow .sh scripts
        d = self.index_input_subreads()
        arrow_scripts = []
        for first_index, last_index in cid_index_chunks:
            arrow_scripts.append(self.create_arrows_bins_no_submit(d, uc, partial_uc,
                                                                 refs, cid_keys[first_index:last_index]))

        # At this point, all arrowed/c{}to{}.sh scripts are ready but not submitted yet
        files_to_run = self.create_submission_jobs(arrow_scripts)
        self.submit_jobs_local_or_remote(files_to_run)

        self.close_log()
        return 0


def add_ice_arrow_arguments(parser):
    """Add arguments for IceArrow2, not including IceArrowPostprocess2."""
    parser = add_cluster_root_dir_as_positional_argument(parser)
    parser = add_fofn_arguments(parser, subread_xml=True)
    parser = add_sge_arguments(parser, arrow_nproc=True, blasr_nproc=True)
    return parser