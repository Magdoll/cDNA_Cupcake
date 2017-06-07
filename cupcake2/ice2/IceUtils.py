"""Define functions useful for IceInit and IceIterative."""
import pdb
import os, sys, subprocess
import os.path as op
import logging
import shutil
import filecmp
import random
import numpy as np
from multiprocessing import Process, Manager
from cPickle import dump, load
from collections import defaultdict
#from pbcore.util.Process import backticks
#from pbcore.io import FastaReader, FastaWriter, FastqReader, FastqWriter, \
#        BasH5Reader
from pbtools.pbtranscript.Utils import realpath, mkdir, \
        get_files_from_fofn, write_files_to_fofn, real_upath
from pbtools.pbtranscript.io.BLASRRecord import BLASRM5Reader
from pbtools.pbtranscript.findECE import findECE
from pbtools.pbtranscript.io.BasQV import basQVcacher
from pbtools.pbtranscript.icedalign.IceDalignUtils import DazzIDHandler, DalignerRunner

__author__ = 'etseng@pacificbiosciences.com'

# define gcon script for ice.
gcon_py = "ice_pbdagcon.py"

# Define data sets for sge sanity check.
dataDir = op.join(op.dirname(op.dirname(op.realpath(__file__))), "data")
GCON_IN_FA = op.join(dataDir, "gcon_in.fa")
GCON_OUT_FA = op.join(dataDir, "gcon_out.fa")

def sanity_check_daligner(scriptDir, testDirName="daligner_test_dir"):
    """
    Run daligner on gcon_in.fa, but don't care about results.
    Just make sure it runs.
    """
    scriptDir = realpath(scriptDir)
    testDir = op.join(scriptDir, testDirName)

    if not op.exists(scriptDir):
        os.makedirs(scriptDir)
    if not op.exists(testDir):
        os.makedirs(testDir)

    testInFa = op.join(testDir, "gcon_in.fa")
    if op.exists(testInFa):
        os.remove(testInFa)
    shutil.copy(GCON_IN_FA, testInFa)
    assert(op.exists(testInFa))

    obj = DazzIDHandler(testInFa)
    DalignerRunner.make_db(obj.dazz_filename)
    runner = DalignerRunner(testInFa, testInFa, is_FL=True, same_strand_only=True, \
                            query_converted=True, db_converted=True, query_made=True, \
                            db_made=True, use_sge=False, cpus=4, sge_opts=None)
    runner.runHPC(min_match_len=300, output_dir=testDir, sensitive_mode=False)

    shutil.rmtree(testDir)
    logging.info("daligner check passed.")
    return True

def sanity_check_gcon():
    """Sanity check gcon."""
    cmd = gcon_py + " --help"
    _out, _code, _msg = backticks(cmd)
    if _code != 0:
        msg = gcon_py + " is not installed."
        raise RuntimeError(msg)
    return gcon_py


def sanity_check_sge(sge_opts, scriptDir, testDirName="gcon_test_dir"):
    """Sanity check if sge can work."""
    scriptDir = realpath(scriptDir)
    testDir = op.join(scriptDir, testDirName)

    if not op.exists(scriptDir):
        os.makedirs(scriptDir)
    if not op.exists(testDir):
        os.makedirs(testDir)

    testSh = op.join(scriptDir, 'test.sh')
    consensusFa = op.join(testDir, "g_consensus.fasta")
    testInFa = op.join(testDir, "gcon_in.fa")
    if op.exists(testInFa):
        os.remove(testInFa)
    shutil.copy(GCON_IN_FA, testInFa)
    assert(op.exists(testInFa))

    with open(testSh, 'w') as f:
        f.write("#!/bin/bash\n")
        f.write("{gcon}".format(gcon=gcon_py) +
                " {inFa} ".format(inFa=real_upath(testInFa)) +
                " {testDir}/g_consensus".format(testDir=real_upath(testDir)) +
                " c1\n")

    assert(op.exists(testSh))
    cmd = "qsub"
    if sge_opts.sge_queue is not None:
        cmd += " -q " + sge_opts.sge_queue
    cmd += " -sync y -pe {env} 1 -cwd -S /bin/bash -V -e /dev/null -o /dev/null {t}".\
          format(t=real_upath(testSh), env=sge_opts.sge_env_name)
    logging.info("Submitting cmd: " + cmd)
    _out, _code, _msg = backticks(cmd)

#    answer = FastaReader(GCON_OUT_FA).__iter__().next()
#    tester = FastaReader(consensusFa).__iter__().next()
#
#    if answer.name != tester.name or \
#       answer.sequence != tester.sequence:
    if not filecmp.cmp(consensusFa, GCON_OUT_FA):
        errMsg = "Trouble running qsub or output is not as " + \
                 "expected ({0} and {1} must agree). Abort!".format(
                     consensusFa, GCON_OUT_FA)
        logging.error(errMsg)
        return False
    else:
        shutil.rmtree(testDir)
        logging.info("sge and gcon check passed.")
        return True



def set_daligner_sensitivity_setting(fastq_filename, is_fasta=False):
    """
    based on the input length range, this sets
     -- sensitivity for DALIGNER  (only two modes: sensitive & NOT sensitive)
     -- how much to ignore in the first 5' end in missed alignment (can override from command line)
    """
    if is_fasta:
        lens = np.array([len(r.sequence) for r in SeqIO.parse(open(fastq_filename), 'fasta')])
    else:
        lens = np.array([len(r.sequence) for r in SeqIO.parse(open(fastq_filename), 'fastq')])
    _low, _high = np.percentile(lens, [1, 95])
    _low  = int(_low)
    _high = int(_high)
    if _low >= 10000:  # for 10kb+
        with open(fastq_filename+'.sensitive.config', 'w') as f:
            f.write("sensitive=True\n")
            f.write("low={0}\nhigh={1}\n".format(_low, _high))
            f.write("ignore5=1000\nignore3=50\n")
            f.write("ece_min_len=60")
            f.write("ece_penalty=1")
        return True, _low, _high, 1000, 50, 60, 1
    elif _low >= 6000:  # for 6-10kb
        with open(fastq_filename+'.sensitive.config', 'w') as f:
            f.write("sensitive=True\n")
            f.write("low={0}\nhigh={1}\n".format(_low, _high))
            f.write("ignore5=800\nignore3=50\n")
            f.write("ece_min_len=40")
            f.write("ece_penalty=1")
        return True, _low, _high, 800, 50, 40, 1
    elif _low >= 3000:  # for 3-6kb
        with open(fastq_filename+'.sensitive.config', 'w') as f:
            f.write("sensitive=False\n")
            f.write("low={0}\nhigh={1}\n".format(_low, _high))
            f.write("ignore5=500\nignore3=50\n")
            f.write("ece_min_len=20")
            f.write("ece_penalty=1")
        return True, _low, _high, 500, 50, 20, 1
    else:# for 0-3 kb
        with open(fastq_filename+'.sensitive.config', 'w') as f:
            f.write("sensitive=False\n")
            f.write("low={0}\nhigh={1}\n".format(_low, _high))
            f.write("ignore5=400\nignore3=50\n")
            f.write("ece_min_len=20")
            f.write("ece_penalty=1")
        return False, _low, _high, 400, 50, 20, 1

def get_daligner_sensitivity_setting(fastq_filename, is_fasta=False):
    """
    config file example (each line must be in order):

    sensitive=False
    low=526
    high=2046
    ignore5=50
    ignore3=50
    ece_min_len=20
    ece_penalty=1
    """
    config = fastq_filename + '.sensitive.config'
    if not os.path.exists(config):
        return set_daligner_sensitivity_setting(fastq_filename, is_fasta)
    else:
        with open(config) as f:
            a, b = f.readline().strip().split('=')
            assert a == 'sensitive'
            flag = (b == 'True')
            a, b = f.readline().strip().split('=')
            assert a == 'low'
            _low = int(b)
            a, b = f.readline().strip().split('=')
            assert a == 'high'
            _high = int(b)
            a, b = f.readline().strip().split('=')
            assert a == 'ignore5'
            _ignore5 = int(b)
            a, b = f.readline().strip().split('=')
            assert a == 'ignore3'
            _ignore3 = int(b)
            a, b = f.readline().strip().split('=')
            assert a == 'ece_min_len'
            _ece_min_len = int(b)
            assert a == 'ece_penalty'
            _ece_penalty = int(b)
        return flag, _low, _high, _ignore5, _ignore3, _ece_min_len, _ece_penalty


def get_ece_arr_from_alignment(record):
    """
    A simplified version of eval_blasr_alignment that does NOT look at QV
    Simply transform record.alnStr to ece_arr where 1 is error

    ex: |||**||||**|||*|*|
    to  000110000110001010
    """
    ece = np.zeros(len(record.alnStr), dtype=np.int)
    for offset, nt_aln in enumerate(record.alnStr):
        ece[offset] = 1 if nt_aln == '*' else 0
    return ece


def eval_blasr_alignment(record, qver_get_func, qvmean_get_func,
                         sID_starts_with_c, qv_prob_threshold, debug=False):
    """
    Takes a BLASRRecord (blasr -m 5) and goes through the
    alignment string
    ex: |||**||||**|||*|*|
    to determine the sequence of 'M' (matches), 'S' (sub), 'I', 'D'

    qver_get_func --- could be either basQV.basQVcacher.get() or
        basQV.basQVcacher.get_smoothed()

    For any non-match, if either or both query/target's QV indicate
    that the event ('S', 'I', 'D') is expected
    (ex: insertion prob >= qv_prob_threshold),
    then it does not count as a penalty.

    Returns: cigar string, binary ECE array

    NOTE: long insertions/deletions are still a difficult problem
    because alignments can be arbitrary right now the quick solution is:
          use probqv get_smoothed
    however, with homopolymers, penalization can still happen unless
    I write code to check specifically for homopolymers, (otherwise the
    cigar_str[-1]=='D' or 'I' sets in). -- Liz
    """
    if debug:
        pdb.set_trace()
    if record.qStrand == '+':
        query_qver_get_func = lambda _name, _pos: qver_get_func(record.qID, _name, min(_pos+record.qStart, record.qLength-1))
        query_qver_get_func1 = lambda _name, _pos: qver_get_func(record.qID, _name, min(_pos+record.qStart+1, record.qLength-1))
    elif record.qStrand == '-':
        query_qver_get_func = lambda _name, _pos: qver_get_func(record.qID, _name, max(0, record.qEnd-1-_pos))
        query_qver_get_func1 = lambda _name, _pos: qver_get_func(record.qID, _name, max(0, record.qEnd-1-_pos-1))
    else:
        raise Exception, "Unknown strand type {0}".format(record.qStrand)

    if record.sStrand == '+':
        subject_qver_get_func = lambda _name, _pos: qver_get_func(record.sID, _name, min(_pos+record.sStart, record.sLength-1))
        subject_qver_get_func1 = lambda _name, _pos: qver_get_func(record.sID, _name, min(_pos+record.sStart+1, record.sLength-1))
    elif record.sStrand == '-':
        subject_qver_get_func = lambda _name, _pos: qver_get_func(record.sID, _name, max(0, record.sEnd-1-_pos))
        subject_qver_get_func1 = lambda _name, _pos: qver_get_func(record.sID, _name, max(0, record.sEnd-1-_pos-1))
    else:
        raise Exception, "Unknown strand type {0}".format(record.sStrand)

    # if mean_qv_for_q|s is not given, always revert back to qv_prob_threshold
    if qvmean_get_func is None:
        mean_qv_for_q = {'D': qv_prob_threshold, 'S': qv_prob_threshold, 'I': qv_prob_threshold}
        mean_qv_for_s = None if sID_starts_with_c else {'D': qv_prob_threshold, 'S': qv_prob_threshold, 'I': qv_prob_threshold}
    else:
        mean_qv_for_q = {'D': qvmean_get_func(record.qID, 'DeletionQV'), \
                         'I': qvmean_get_func(record.qID, 'InsertionQV'), \
                         'S': qvmean_get_func(record.qID, 'SubstitutionQV')}
        if sID_starts_with_c:
            mean_qv_for_s = None
        else:
            mean_qv_for_s = {'D': qvmean_get_func(record.sID, 'DeletionQV'), \
                            'I': qvmean_get_func(record.sID, 'InsertionQV'), \
                            'S': qvmean_get_func(record.sID, 'SubstitutionQV')}

    q_index = 0
    s_index = 0
    last_state, last_tracking_nt, homopolymer_so_far = None, None, False
    cigar_str = ''
    # binary array of 0|1 where 1 is a penalty
    ece = np.zeros(len(record.alnStr), dtype=np.int)
    # pdb.set_trace()
    for offset, nt_aln in enumerate(record.alnStr):
        if nt_aln == '|':  # match
            cigar_str += 'M'
            q_index += 1
            s_index += 1
            last_state = 'M'
        elif record.qAln[offset] == '-':  # deletion
            # if last_state != 'D':
            #     last_tracking_nt = record.sAln[s_index]
            #     homopolymer_so_far = True
            #     ece[offset] = 1
            # else:
            #     homopolymer_so_far = (record.sAln[s_index] == last_tracking_nt)
            #     if not homopolymer_so_far:
            #         ece[offset] = 1

            # for deletion, cases where consider a non-match
            # (1) IF last position was not "D"
            #        case 1a: both query and subject has very good prob
            # (2) IF last position was "D" (q_index did not advance)
            #        case 2a: both query and subject has very good prob
            #        case 2b: subject has good prob; query has bad prob AND last-cur pos is NOT homopolymer
            # case 1a and case 2a have the same condition
            # case 2b arises because the only possible explanation would have been query have bad prob,
            #   but it was used to explain the last deletion and the advanced S nucleotide is diff from the last S
            s_is_good = sID_starts_with_c or subject_qver_get_func('InsertionQV', s_index) < mean_qv_for_s
            q_is_good = query_qver_get_func1('DeletionQV', q_index) < mean_qv_for_q
            if last_state != 'D': # entering D state now, record s
                last_tracking_nt = record.sAln[offset]
                homopolymer_so_far = True
                if (s_is_good and q_is_good):
                    ece[offset] = 1
            else:  # already in D state, which means q_index did not advance, hence q_is_good is the same
                homopolymer_so_far = (record.sAln[offset] == last_tracking_nt)
                if (s_is_good and (q_is_good or not homopolymer_so_far)):
                    ece[offset] = 1


            # ------------- OLD VERSION -------------------
            #if (offset > 0 and cigar_str[-1] == 'D') or (
            #    query_qver_get_func1('DeletionQV', q_index) < qv_prob_threshold and
            #   (sID_starts_with_c or subject_qver_get_func('InsertionQV', s_index) < qv_prob_threshold)):
            #    # case 1: last one was also a D (so q_index did not advance)
            #    # case 2: both QVs were good yet still a non-match, penalty!
            #    ece[offset] = 1
            # ------------- OLD VERSION -------------------
            cigar_str += 'D'
            s_index += 1
            last_state = 'D'
        elif record.sAln[offset] == '-':  # insertion
            # if last_state != 'I':
            #     last_tracking_nt = record.qAln[q_index]
            #     homopolymer_so_far = True
            #     ece[offset] = 1
            # else:
            #     homopolymer_so_far = (record.qAln[q_index] == last_tracking_nt)
            #     if not homopolymer_so_far:
            #         ece[offset] = 1

            # for insertion, cases where consider a non-match
            # (1) IF last position was not "I"
            #        case 1a: both query and subject has very good prob
            # (2) IF last position was "I" (s_index did not advance)
            #        case 2a: both query and subject has very good prob
            #        case 2b: query has good prob; subject has bad prob AND last-cur pos is NOT homopolymer
            q_is_good = query_qver_get_func('InsertionQV', q_index) < mean_qv_for_q
            s_is_good = sID_starts_with_c or subject_qver_get_func1('DeletionQV', s_index) < mean_qv_for_s

            if last_state != 'I':
                last_tracking_nt = record.qAln[offset]
                homopolymer_so_far = True
                if (q_is_good and s_is_good):
                    ece[offset] = 1
            else: # already in "I" state, s_index did not advance, s_is_good is the same
                homopolymer_so_far = (record.qAln[offset] == last_tracking_nt)
                if q_is_good and (s_is_good or not homopolymer_so_far):
                    ece[offset] = 1

            # ------------- OLD VERSION -------------------
            #if (offset > 0 and cigar_str[-1] == 'I') or (
            #        query_qver_get_func('InsertionQV', q_index) < qv_prob_threshold and
            #        (sID_starts_with_c or subject_qver_get_func1('DeletionQV', s_index) < qv_prob_threshold)):
            #    # case 1: last one was also a I (so s_index did not advance)
            #    # case 2: both QVs were good yet still a no-match
            #    ece[offset] = 1
            # ------------- OLD VERSION -------------------
            cigar_str += 'I'
            q_index += 1
            last_state = 'I'
        else:  # substitution
            cigar_str += 'S'
            if query_qver_get_func('SubstitutionQV', q_index) < mean_qv_for_q and \
               (sID_starts_with_c or subject_qver_get_func('SubstitutionQV', s_index) < mean_qv_for_s):
                ece[offset] = 1
            q_index += 1
            s_index += 1
            last_state = 'S'


    return cigar_str, ece


class HitItem(object):

    """
    Simply define an object class for saving items produced by
    blasr_against_ref.
    """

    def __init__(self, qID, cID, qStart=None, qEnd=None,
                 missed_q=None, missed_t=None,
                 fakecigar=None, ece_arr=None):
        self.qID = qID
        self.cID = cID
        self.qStart = qStart
        self.qEnd = qEnd
        self.missed_q = missed_q
        self.missed_t = missed_t
        self.fakecigar = fakecigar
        self.ece_arr = ece_arr


def blasr_against_ref(output_filename, is_FL, sID_starts_with_c,
                      qver_get_func, qvmean_get_func, qv_prob_threshold=.03,
                      ece_penalty=1, ece_min_len=20, same_strand_only=True, max_missed_start=50, max_missed_end=50):
    """
    Excluding criteria:
    (1) self hit
    (2) opposite strand hit  (should already be in the same orientation;
        can override with <same_strand_only> set to False)
    (3) less than 90% aligned or more than 50 bp missed

    qver_get_func --- should be basQV.basQVcacher.get() or
                      .get_smoothed(), or can just pass in
                      lambda (x, y): 1. to ignore QV
    """
    with BLASRM5Reader(output_filename) as reader:
        for r in reader:
            missed_q = r.qStart + r.qLength - r.qEnd
            missed_t = r.sStart + r.sLength - r.sEnd

            if sID_starts_with_c:
                # because all consensus should start with
                # c<cluster_index>
                assert r.sID.startswith('c')
                if r.sID.find('/') > 0:
                    r.sID = r.sID.split('/')[0]
                if r.sID.endswith('_ref'):
                    # probably c<cid>_ref
                    cID = int(r.sID[1:-4])
                else:
                    cID = int(r.sID[1:])
            else:
                cID = r.sID

            # self hit, useless!
            # low identity not allowed
            # opposite strand not allowed!
            if (cID == r.qID or
                    r.identity < 70. or
                    (r.strand == '-' and same_strand_only)):
                yield HitItem(qID=r.qID, cID=cID)
                continue

            # full-length case: allow up to <max_missed_start> bp of 5' not aligned
            # and 50bp of 3' not aligned
            # non-full-length case: not really tested...don't use
            if is_FL and (r.sStart > max_missed_start or r.qStart > max_missed_start or
                          (r.sLength - r.sEnd > max_missed_end) or
                          (r.qLength - r.qEnd > max_missed_end)):
                yield HitItem(qID=r.qID, cID=cID)
            else:
                cigar_str, ece_arr = eval_blasr_alignment(
                    record=r,
                    qver_get_func=qver_get_func,
                    sID_starts_with_c=sID_starts_with_c,
                    qv_prob_threshold=qv_prob_threshold,
                    qvmean_get_func=qvmean_get_func)

                if alignment_has_large_nonmatch(ece_arr,
                                                ece_penalty, ece_min_len):
                    yield HitItem(qID=r.qID, cID=cID)
                else:
                    yield HitItem(qID=r.qID, cID=cID,
                                  qStart=r.qStart, qEnd=r.qEnd,
                                  missed_q=missed_q * 1. / r.qLength,
                                  missed_t=missed_t * 1. / r.sLength,
                                  fakecigar=cigar_str,
                                  ece_arr=ece_arr)


def alignment_has_large_nonmatch(ece_arr, penalty, min_len):
    """
    penalty of (-)1: 50%
    penalty of (-)2: 66%
    penalty of (-)4: 80%
    penalty of (-)9: 90%

    Return True when alignment has large non-matches not explained
    by low base QVs (in other words, "reject" as an isoform hit and
    don't put in the same cluster)
    """
    ece_arr = ece_arr * (penalty + 1)
    s = [0] + list(ece_arr - penalty)
    # fix this later to something faster & better
    return (len(findECE(s, len(s), min_len, True)) > 0)


def possible_merge(r, ece_penalty, ece_min_len, max_missed_start, max_missed_end):
    """
    r --- BLASRM5Record
    Criteria:
    (1) identity >= 90% and same strand
    (2) check criteria for how much is allowed to differ on the
        5' / 3' ends
    """
    if r.sID == r.qID or r.identity < 90 or r.strand == '-':
        return False
    # intentional here to prevent disrupting future ICE runs
    # MORE lenient on 5' but NOT on 3'
    if ((r.qLength - r.qEnd) > max_missed_end or (r.sLength - r.sEnd) > max_missed_end or
            r.qStart > max_missed_start or r.sStart > max_missed_start):
        return False

    arr = np.array([(x == '*') * 1 for x in r.alnStr])
    if alignment_has_large_nonmatch(ece_arr=arr,
                                    penalty=ece_penalty,
                                    min_len=ece_min_len):
        return False
    return True


def get_the_only_fasta_record(fa):
    """Input fasta file should contain exactly one FastaRecord,
    return the fastas record."""
    rs = [r for r in FastaReader(fa)]
    if len(rs) != 1:
        errMsg = "Cluster fasta file {fa} must contain only one read.".\
            format(fa=fa)
        raise ValueError(errMsg)
    return rs[0]


"""
The following methods was originally created by jchin:
    ~jchin/depot_mp27/jchin/rset_quvier.py
, and then modified by etseng.

Input: input.fasta.fofn (shared),
       per-cluster in.fa,
       per-cluster g_consensus.fa

-- input.fasta.fofn should be raw fasta files
   (pls2fasta -maskRegion) of input.fofn

Within each cluster:
1) create in.raw.fa based on input.fasta.fofn & in.fa,
   putting in raw (unrolled) fasta of each ZMW
2) blasr (1) to g_consensus.fa output as SAM

This is faster than using regions.fofn because it still reads
through the whole .bax.h5 files
"""


def is_blank_sam(samfile):
    """
    return True if the SAM file only has @xx header and NO alignment
    """
    with open(samfile) as f:
        for line in f:
            if not line.startswith('@'):
                return False
    return True


def concat_sam(samfiles, outsam_filename):
    """
    Header looks like:
    @HD     VN:1.3.1
    @SQ     SN:c31  LN:3104 M5:ef7d3f84dea9d9face43e6fd5b6336c4
    @RG     ID:2caa54eef6   PU:in.raw_with_partial.fa       SM:NO_CHIP_ID
    @PG     ID:BLASR        VN:1.3.1.126469 CL:blasr in.raw_with_partial.fa g_consensus.fa -nproc 12 -bestn 5 -nCandidates 10 -sam -out out.sam

    NOTE: check for M5 conflicts; manipulate them if it conflicts
    """
    f_sq = open(outsam_filename + '.sq', 'w')
    f_bd = open(outsam_filename + '.bd', 'w')

    rg_line = None
    pg_line = None

    md5_seen = set()

    if len(samfiles) == 0:
        raise ValueError("No sam input files to concatenate.")

    h = open(samfiles[0])
    line = h.readline()
    assert line.startswith('@HD')
    f_sq.write(line)
    line = h.readline()
    assert line.startswith('@SQ')
    line = h.readline()
    assert line.startswith('@RG')
    rg_line = line  # write at the end
    line = h.readline()
    assert line.startswith('@PG')
    pg_line = line  # write at the end
    h.close()

    for f in samfiles:
        with open(f) as h:
            assert h.readline().startswith('@HD')
            line = h.readline()
            assert line.startswith('@SQ')
            # ------- check for MD5 conflicts ----------- #
            m5 = line.strip().split()[-1]
            assert m5.startswith("M5:")
            if m5 not in md5_seen:
                f_sq.write(line)
                md5_seen.add(m5)
            else:
                s = list(m5[3:])
                while True:
                    # create a random m5 string.
                    random.shuffle(s)
                    s = "".join(s)
                    if s not in md5_seen:
                        break
                line = line[:line.find('M5:')] + 'M5:' + s + '\n'
                logging.debug("MD5 conflict: change to {0}".format(s))
                md5_seen.add(s)
                f_sq.write(line)
            # ----- end MD5 checking and writing --------- #
            assert h.readline().startswith('@RG')
            assert h.readline().startswith('@PG')
            for line in h:
                f_bd.write(line)

    f_bd.close()
    f_sq.write(rg_line)
    f_sq.write(pg_line)
    f_sq.close()

    cmd = "cat {0}.sq {0}.bd > {0}".format(real_upath(outsam_filename))
    _out, _code, _msg = backticks(cmd)
    if _code != 0:
        raise IOError("Failed to concat sam files! Abort." + _msg)

    os.remove(f_sq.name)
    os.remove(f_bd.name)

def convert_fofn_to_fasta_worker(in_queue):
    while not in_queue.empty():
        stuff = in_queue.get()
        pls2fasta_cmd, tmp_out_file, out_file = stuff
        _out, _code, _msg = backticks(pls2fasta_cmd)
        if _code != 0:
            raise RuntimeError("CMD failed: {cmd}\n".format(cmd=pls2fasta_cmd) + _msg)
        trim_subread_flanks(tmp_out_file, out_file)
        if op.exists(tmp_out_file):
            os.remove(tmp_out_file)

def convert_fofn_to_fasta(fofn_filename, out_filename, fasta_out_dir,
                          force_overwrite=False, cpus=1):
    """
    For each .bax.h5 file, create .bax.h5.fasta file and save paths to
    out_filename, which should usually be 'input.fasta.fofn'
    """
    logging.info("Converting fofn {fofn} to fasta.".format(fofn=fofn_filename))
    in_fns = get_files_from_fofn(fofn_filename)
    #out_fns = []
    mkdir(fasta_out_dir)

    # multiprocessing worker stuff
    manager = Manager()
    in_queue = manager.Queue(len(in_fns))
    in_queue_count = 0
    outfile_track = {} # expected out file --> (cmd, tmp)
    pool = []
    out_fns = []

    for in_fn in in_fns:
        #print >> sys.stderr, "DEBUG: converting h5 file:", in_fn
        logging.debug("converting h5 file: {f}.".format(f=in_fn))
        if not (in_fn.endswith('.bax.h5') or in_fn.endswith('.bas.h5')):
            raise ValueError("fofn file {fofn} ".format(fofn=fofn_filename) +
                             "should only contain bax/bas.h5 files.")

        # e.g. m111xxxx.1.bax.h5 ==>
        #      tmp_out_file = m11xxxx.1.bax.h5.fasta.tmp
        #      out_file = m11xxxx.1.bax.h5.fasta
        in_basename = op.basename(in_fn)
        tmp_out_file = op.join(fasta_out_dir, in_basename + '.fasta.tmp')
        out_file = op.join(fasta_out_dir, in_basename + '.fasta')
        if op.exists(out_file) and not force_overwrite:
            logging.debug("File {0} already exists. skipping.".format(out_file))
            out_fns.append(out_file)
            if op.exists(tmp_out_file):
                os.remove(tmp_out_file)
        else:
            cmd = "pls2fasta {in_fn} ".format(in_fn=real_upath(in_fn)) + \
                  " {out} ".format(out=real_upath(tmp_out_file)) + \
                  "-minSubreadLength 300 -minReadScore 750 -trimByRegion"
            print >> sys.stderr, "DEBUG: putting in queue:", cmd, tmp_out_file, out_file
            in_queue.put((cmd, tmp_out_file, out_file))
            in_queue_count += 1
            outfile_track[out_file] = (cmd, tmp_out_file)
            print >> sys.stderr, "DEBUG: put in queue:", cmd, tmp_out_file, out_file

    cpus = min(cpus, in_queue_count) # cap max CPU if there's fewer files to convert
    for i in xrange(cpus):
        p = Process(target=convert_fofn_to_fasta_worker, args=(in_queue,))
        pool.append(p)

    #error_flag = False
    # starting & joining pool worakers
    for p in pool:
        p.start()
        #print >> sys.stderr, "Starting worker", p.name
    for p in pool:
        #print >> sys.stderr, "Waiting join", p.name
        p.join(timeout=1200)
        if p.is_alive(): p.terminate()

    # check that all files exists
    # if it does not, force to run locally
    for out_file,(cmd, tmp_out_file) in outfile_track.iteritems():
        in_queue.put((cmd, tmp_out_file, out_file))
        convert_fofn_to_fasta_worker(in_queue)
        out_fns.append(out_file)

    #if error_flag:
    #    raise Exception, "Unable to successfuly run convert_fofn_to_fasta, ABORT!"

    write_files_to_fofn(out_fns, out_filename)


def trim_subread_flanks(fasta_filename, output_filename,
                        trim_len=100, min_len=100):
    """
    fasta_filename --- should be subread output from pls2fasta

    trim first/last 100bp (which contains primer&polyA) away and correct
    coordinates
    """
    with FastaWriter(output_filename) as writer, \
            FastaReader(fasta_filename) as reader:
        for r in reader:
            # ex: m14011..._s1_p0/15/1305_4354
            movie, hn, s_e = r.name.split()[0].split('/')
            s, e = s_e.split('_')
            s, e = int(s), int(e)
            assert s < e
            s2 = s + trim_len
            e2 = e - trim_len
            if e2 - s2 >= min_len:
                newname = "{0}/{1}/{2}_{3}".format(movie, hn, s2, e2)
                newseq = r.sequence[trim_len:-trim_len]
                writer.writeRecord(newname, newseq)


def build_sa(input_fasta, out_sa):
    """Generate suffix array of input_fasta"""
    if op.exists(input_fasta):
        cmd = "sawriter {o} {i} -blt 8 -welter ".\
            format(o=real_upath(out_sa), i=real_upath(input_fasta))
        _out, _code, _msg = backticks(cmd)
        if _code == 0:
            return True
        else:
            # If failed to generate suffix array, warning.
            logging.warn("Unable to create suffix array for {f}.".format(f=input_fasta))
            return False
    else:
        raise IOError("Unable to find fasta file {f}.".format(f=input_fasta))


def write_in_raw_fasta_starhelper(args):
    write_in_raw_fasta(*args)

def write_in_raw_fasta(input_fasta_d, in_seqids,
                       out_fa, ignore_keyerror):
    """
    input_fasta_d --- miscBio.MetaFastaReader
    input fasta should be in format <movie>/<holeNumber>/<subread or CCS stuff>

    Create a out_fa fasta where we dump the "raw" (unrolled) of every ZMW from in_fa
    """
    movies = set()
    zmw_seen = set()
    with open(out_fa, 'w') as f:
        for seqid in in_seqids:
            try:
                zmw = seqid[:seqid.rfind('/')]
                if zmw not in zmw_seen:
                    movies.add(zmw.split('/')[0])
                    for rec in input_fasta_d[zmw]:
                        f.write(">{0}\n{1}\n".format(rec.name, rec.sequence))
                    zmw_seen.add(zmw)
            except KeyError:
                if ignore_keyerror:
                    pass#logging.warning("Ignoring {zmw} because the input fasta_fofn ".
                        #            format(zmw=zmw) + "does not contain its sequence.")
                else:
                    raise ValueError, "{0} doesn't exist. Abort!".format(zmw)
    #return movies


def blasr_sam_for_quiver(input_fasta, ref_fasta,
                         out_sam_filename,
                         run_cmd=True, blasr_nproc=12):
    """
    input_fasta --- should be in.raw.fa
    ref_fasta --- reference fasta (ex: g_consensus.fa) to align to
    out_sam_filename --- sam output aligning in_fasta to ref_fasta

    run blasr -clipping soft to get sam
    """
    cmd = "blasr {i} ".format(i=real_upath(input_fasta)) + \
          "{r} ".format(r=real_upath(ref_fasta)) + \
          "-nproc {n} ".format(n=blasr_nproc) + \
          "-bestn 5 -nCandidates 10 -sam -clipping soft " + \
          "-out {o}".format(o=real_upath(out_sam_filename))
    logging.debug("CMD: " + cmd)
    if run_cmd:
        _out, _code, _msg = backticks(cmd)
        if _code != 0:
            raise RuntimeError("CMD failed: {cmd}\n{e}".
                               format(cmd=cmd, e=_msg))
    return cmd


def num_reads_in_fasta(in_fa):
    """Return the number of reads in the in_fa fasta file."""
    if not op.exists(in_fa):
        raise IOError("fasta file {f} does not exist.".format(f=in_fa))
    cmd = "grep '>' {f} | wc -l ".format(f=in_fa)
    logging.debug("CMD: " + cmd)
    _out, _code, _msg = backticks(cmd)
    if _code != 0:
        raise RuntimeError("CMD failed: {cmd}\n{e}".
                           format(cmd=cmd, e=_msg))
    return int(_out[0])


def combine_nfl_pickles(splitted_pickles, out_pickle):
    """Combine splitted nfl pickles to a big pickle."""
    logging.debug("Cominbing {N} nfl pickles: {ps} ".
                  format(N=len(splitted_pickles),
                         ps=",".join(splitted_pickles)) +
                  " into a big pickle {p}.".format(p=out_pickle))

    if len(splitted_pickles) == 1:
        logging.debug("Copying the only given pickle to out_pickle.")
        if realpath(splitted_pickles[0]) != realpath(out_pickle):
            shutil.copyfile(splitted_pickles[0], out_pickle)
    else:
        # Combine all partial outputs
        logging.debug("Merging all pickles.")
        partial_uc = defaultdict(lambda: [])
        nohit = set()
        for pf in splitted_pickles:
            logging.debug("Merging {pf}.".format(pf=pf))
            a = load(open(pf))
            nohit.update(a['nohit'])
            for k, v in a['partial_uc'].iteritems():
                partial_uc[k] += v

        logging.debug("Dumping all to {f}".format(f=out_pickle))
        # Dump to one file
        partial_uc = dict(partial_uc)
        with open(out_pickle, 'w') as f:
            dump({'nohit': nohit, 'partial_uc': partial_uc}, f)
        logging.debug("{f} created.".format(f=out_pickle))

def cid_with_annotation(cid):
    """Given a cluster id, return cluster id with human readable annotation.
    e.g., c0 --> c0 isoform=c0
          c0/89/3888 -> c0/89/3888 isoform=c0;full_length_coverage=89;isoform_length=3888
          c0/f89p190/3888 -> c0/f89p190/3888 isoform=c0;full_length_coverage=89;non_full_length_coverage=190;isoform_length=3888
    """
    fields = cid.split('/')
    short_id, fl_coverage, nfl_coverage, seq_len = None, None, None, None
    if len(fields) != 1 and len(fields) != 3:
        raise ValueError("Not able to process isoform id: {cid}".format(cid=cid))
    short_id = fields[0]
    if len(fields) == 3:
        seq_len = fields[2]
        if "f" in fields[1]:
            if "p" in fields[1]: # f89p190
                fl_coverage = fields[1].split('p')[0][1:]
                nfl_coverage = fields[1].split('p')[1]
            else: # f89
                fl_coverage = fields[1][1:]
        else:
            fl_coverage = fields[1]

    annotations = ["isoform={short_id}".format(short_id=short_id)]
    if fl_coverage is not None:
        annotations.append("full_length_coverage={fl}".format(fl=fl_coverage))
    if nfl_coverage is not None:
        annotations.append("non_full_length_coverage={nfl}".format(nfl=nfl_coverage))
    if seq_len is not None:
        annotations.append("isoform_length={l}".format(l=seq_len))

    return "{cid} {annotation}".format(cid=cid, annotation=";".join(annotations))


def get_qv_from_bas_handler(bas_handler, hn, s_e, qv_name):
    """Read QV of type qv_name for movie/hn/s_e from input bas_h5 file handler."""
    is_CCS = True  # Assume read is CCS
    strand = '+'
    s, e = 0, 0
    get_all = False
    try:
        if s_e == "ccs":
            is_CCS = True
            get_all = True
        elif s_e.endswith('_CCS'):
            is_CCS = True
            s, e = s_e.split('_')[:2]
        else:
            is_CCS = False
            s, e = s_e.split('_')
    except ValueError:
        raise ValueError("{s_e} is not a valid read start_end.".
                         format(s_e=s_e))
    s, e = int(s), int(e)
    if s > e:
        s, e = e, s
        strand = '-'

    zmw = bas_handler[hn]
    if is_CCS:
        if zmw.ccsRead is not None:
            qvs = zmw.ccsRead.qv(qv_name)
        else:  # this is for reads_of_insert w/ 0-passed
            qvs = zmw.read().qv(qv_name)
    else:  # subread
        qvs = zmw.read().qv(qv_name)
    if not get_all:
        qvs = qvs[s:e]
    if strand == '-':
        qvs = qvs[::-1]
    return qvs

def ice_fq2fa(in_fq, out_fa):
    handle = FastaWriter(out_fa)
    for r in FastqReader(in_fq):
        handle.writeRecord(r.name, r.sequence)

def ice_fa2fq(in_fa, ccs_fofn, out_fq):
    """Convert an input FASTA file to an output FASTQ file,
       reading QVs from the input ccs.h5 or ccs FOFN.
    """

    qver = basQVcacher()
    if ccs_fofn.endswith(".h5"):  # Input is a ccs.h5 file not a FOFN.
        qver.add_bash5(ccs_fofn)
    else:  # Input is a ccs FOFN containing multiple ccs.h5 files.
        for ccs_fn in get_files_from_fofn(ccs_fofn):
            qver.add_bash5(ccs_fn)

    bas_handlers = {}

    with FastaReader(in_fa) as reader, \
            FastqWriter(out_fq) as writer:
        for r in reader:
            seqid = r.name.split(' ')[0]
            movie, hn, s_e = "", "", ""
            try:
                movie, hn, s_e = seqid.split('/')
                hn = int(hn)
            except ValueError:
                raise ValueError("{seqid} is not a valid CCS read".
                                 format(seqid=seqid))
            try:
                bas_file = qver.bas_files[movie][seqid]
                if bas_file not in bas_handlers:
                    bas_handlers[bas_file] = BasH5Reader(bas_file)
            except KeyError:
                raise IOError("Could not read {s} from input ccs fofn.".
                              format(s=seqid))
            logging.debug("Getting QVs for {name} ...".format(name=r.name))
            qvs = get_qv_from_bas_handler(bas_handler=bas_handlers[bas_file],
                                          hn=hn, s_e=s_e,
                                          qv_name="QualityValue")
            if len(r.sequence) != len(qvs):
                raise ValueError("Sequence and QVs of {r} should be the same!".
                                 format(r=r.name))
            writer.writeRecord(r.name, r.sequence, qvs)

    for bas_file, bas_handler in bas_handlers.iteritems():
        logging.debug("Closing {bas_file} ...".format(bas_file=bas_file))
        bas_handler.close()


def locally_run_failed_quiver_jobs(bad_sh_files, max_fail=3):
    """
    bad_sh_files --- list of xxx.sh Quiver jobs that failed to run. Attempt to re-run.
    """
    still_bad = []
    for sh_file in bad_sh_files:
        err_count = 0
        failed = True
        while err_count <= max_fail:
            try:
                subprocess.check_call("bash " + sh_file, shell=True)
                failed = True
                break
            except:
                err_count += 1
        if failed: still_bad.append(sh_file)
    return still_bad

