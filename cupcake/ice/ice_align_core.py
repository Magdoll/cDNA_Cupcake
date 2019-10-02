__author__ = 'etseng@pacb.com'

"""
Replicate of key methods from ICE related to alignment evaluation and ECE
"""

import pdb
import numpy as np
from cupcake.ice.find_ECE import findECE

def get_ece_arr_from_alignment(record):
    """
    A simplified version of eval_sam_alignment that does NOT look at QV
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
        raise Exception("Unknown strand type {0}".format(record.qStrand))

    if record.sStrand == '+':
        subject_qver_get_func = lambda _name, _pos: qver_get_func(record.sID, _name, min(_pos+record.sStart, record.sLength-1))
        subject_qver_get_func1 = lambda _name, _pos: qver_get_func(record.sID, _name, min(_pos+record.sStart+1, record.sLength-1))
    elif record.sStrand == '-':
        subject_qver_get_func = lambda _name, _pos: qver_get_func(record.sID, _name, max(0, record.sEnd-1-_pos))
        subject_qver_get_func1 = lambda _name, _pos: qver_get_func(record.sID, _name, max(0, record.sEnd-1-_pos-1))
    else:
        raise Exception("Unknown strand type {0}".format(record.sStrand))

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