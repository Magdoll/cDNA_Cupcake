
import os, re, sys, subprocess
import numpy as np
from cupcake.io.BioReaders import GMAPSAMReader
from cupcake.ice.find_ECE import findECE

gcon2_py = "ice_pbdagcon2.py"


def sanity_check_gcon2():
    """Sanity check gcon."""
    cmd = gcon2_py + " --help"

    errmsg = gcon2_py + " is not installed."
    if subprocess.check_call(cmd, shell=True)!=0:
        print >> sys.stderr, "ERROR RUNNING CMD:", cmd
        print >> sys.stderr, errmsg
        sys.exit(-1)
    return gcon2_py


def alignment_missed_start_end_less_than_threshold(r, max_missed_start, max_missed_end,
                   full_missed_start, full_missed_end):
    """
    Check that whichever is the shorter one, must be close to fully mapped
    (subject to full_missed_start/end)
    and the longer one is allowed to have more missed start/end
    (subject to max_missed_start/end)

    """
    assert max_missed_start >= full_missed_start and max_missed_end >= full_missed_end

    # which ever is the shorter one, must be fully mapped
    missed_start_1 = r.qStart
    missed_start_2 = r.sStart
    missed_end_1 = (r.qLen - r.qEnd)
    missed_end_2 = (r.sLen - r.sEnd)
    if r.qLen > r.sLen:
        missed_start_1, missed_start_2 = missed_start_2, missed_start_1
        missed_end_1, missed_end_2 = missed_end_2, missed_end_1
    # the smaller one must be close to fully mapped
    if (missed_start_1 > full_missed_start) or \
            (missed_end_1 > full_missed_end) or \
            (missed_start_2 > max_missed_start) or \
            (missed_end_2 > max_missed_end):
        return False

    return True

def minimap2_against_ref2(sam_filename, query_len_dict, ref_len_dict,
                      is_FL, sID_starts_with_c,
                      ece_penalty=1, ece_min_len=20, same_strand_only=True,
                      max_missed_start=200, max_missed_end=50,
                      full_missed_start=50, full_missed_end=30):
    """
    Excluding criteria:
    (1) self hit
    (2) opposite strand hit  (should already be in the same orientation;
        can override with <same_strand_only> set to False)
    """
    for r in GMAPSAMReader(sam_filename, True, query_len_dict=query_len_dict, ref_len_dict=ref_len_dict):
            missed_q = r.qStart + r.qLen - r.qEnd
            missed_t = r.sStart + r.sLen - r.sEnd

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
            # opposite strand not allowed!
            if (cID == r.qID or (r.flag.strand == '-' and same_strand_only)):
                yield HitItem(qID=r.qID, cID=cID)
                continue

            # regardless if whether is full-length (is_FL)
            # the query MUST be mapped fully (based on full_missed_start/end)
            if r.qStart > full_missed_start or (r.qLen-r.qEnd) > full_missed_end:
                yield HitItem(qID=r.qID, cID=cID)

            # full-length case: allow up to max_missed_start bp of 5' not aligned
            # and max_missed_end bp of 3' not aligned
            # non-full-length case: not really tested...don't use
            if is_FL and not alignment_missed_start_end_less_than_threshold(r,\
                            max_missed_start, max_missed_end, full_missed_start, full_missed_end):
                yield HitItem(qID=r.qID, cID=cID)
            else:
                ece_arr = eval_sam_alignment(r)

                if alignment_has_large_nonmatch(ece_arr,
                                                ece_penalty, ece_min_len):
                    yield HitItem(qID=r.qID, cID=cID)
                else:
                    yield HitItem(qID=r.qID, cID=cID,
                                  qStart=r.qStart, qEnd=r.qEnd,
                                  missed_q=missed_q * 1. / r.qLen,
                                  missed_t=missed_t * 1. / r.sLen,
                                  fakecigar=r.cigar,
                                  ece_arr=ece_arr)


def possible_merge2(r, ece_penalty, ece_min_len,
                   max_missed_start=200, max_missed_end=50,
                   full_missed_start=50, full_missed_end=30):
    """
    r --- BLASRM5Record
    Criteria:
    (1) identity >= 90% and same strand
    (2) check criteria for how much is allowed to differ on the
        5' / 3' ends

    Note: one must be fully mapped (allowing only a small portion to be unmapped)
          while the other can have <max_missed_start>/<max_missed_end>
    """
    if r.sID == r.qID or r.identity < 90 or r.strand == '-':
        return False


    if not alignment_missed_start_end_less_than_threshold(r, max_missed_start, max_missed_end,
                                            full_missed_start, full_missed_end):
        return False

    arr = np.array([(x == '*') * 1 for x in r.alnStr])
    if alignment_has_large_nonmatch(ece_arr=arr,
                                    penalty=ece_penalty,
                                    min_len=ece_min_len):
        return False
    return True


def cid_with_annotation2(cid, expected_acc=None):
    """Given a cluster id, return cluster id with human readable annotation.
    e.g., c0 --> c0 isoform=c0
          c0/89/3888 -> c0/89/3888 isoform=c0;full_length_coverage=89;isoform_length=3888;expected_accuracy=0.99
          c0/f89p190/3888 -> c0/f89p190/3888 isoform=c0;full_length_coverage=89;non_full_length_coverage=190;isoform_length=3888;expected_accuracy=0.99
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
    if expected_acc is not None:
        annotations.append("expected_accuracy={0:.3f}".format(expected_acc))

    return "{cid} {annotation}".format(cid=cid, annotation=";".join(annotations))


def eval_sam_alignment(record, debug=False):
    """
    Takes a GMAPSAMRecord and goes through the cigar string

    ToDo: take in QV information & account homopolymer in the future

    Returns: binary ECE array
    """
    if debug:
        import pdb
        pdb.set_trace()

    # binary array of 0|1 where 1 is a penalty
    aln_len = record.num_mat_or_sub + record.num_ins + record.num_del
    ece = np.zeros(aln_len, dtype=np.int)

    q_index = record.qStart
    s_index = record.sStart
    offset = 0
    for _match in re.finditer('(\d+)(\S)', record.cigar):
        _count, _type = _match.groups()
        _count = int(_count)
        if _type in ('M', '='):
            q_index += _count
            s_index += _count
            offset += _count
        elif _type == 'X':  # mismatch
            ece[offset:offset+_count] = 1
            q_index += _count
            s_index += _count
            offset += _count
        elif _type == 'D':
            ece[offset:offset+_count] = 1
            s_index += _count
            offset += _count
        elif _type == 'I':
            ece[offset:offset+_count] = 1
            q_index += _count
            offset += _count
        elif _type in ('H', 'S'):
            q_index += _count
        else:
            raise Exception, "Unrecognized cigar character {0}!".format(_type)


    return ece

class HitItem(object):

    """
    Simply define an object class for saving items produced by
    blasr_against_ref or daligner_against_ref.
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

    def __str__(self):
        return """{qID}/{qStart}_{qEnd} aligns to {cID}""".format(
                qID=self.qID.split(' ')[0], cID=self.cID.split(' ')[0],
                qStart=self.qStart, qEnd=self.qEnd)

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
    return len(findECE(s, len(s), min_len, True)) > 0