
import numpy as np
from pbtranscript.Utils import execute
from pbtranscript.ice.IceUtils import alignment_has_large_nonmatch, HitItem, eval_blasr_alignment
from pbtranscript.io import LA4IceReader, BLASRM5Reader

gcon2_py = "ice_pbdagcon2.py"

def sanity_check_gcon2():
    """Sanity check gcon."""
    cmd = gcon2_py + " --help"

    errmsg = gcon2_py + " is not installed."
    execute(cmd=cmd, errmsg=errmsg)
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
    missed_end_1 = (r.qLength - r.qEnd)
    missed_end_2 = (r.sLength - r.sEnd)
    if r.qLength > r.sLength:
        missed_start_1, missed_start_2 = missed_start_2, missed_start_1
        missed_end_1, missed_end_2 = missed_end_2, missed_end_1
    # the smaller one must be close to fully mapped
    if (missed_start_1 > full_missed_start) or \
            (missed_end_1 > full_missed_end) or \
            (missed_start_2 > max_missed_start) or \
            (missed_end_2 > max_missed_end):
        return False

    return True

def blasr_against_ref2(output_filename, is_FL, sID_starts_with_c,
                      qver_get_func, qvmean_get_func, qv_prob_threshold=.03,
                      ece_penalty=1, ece_min_len=20, same_strand_only=True,
                      max_missed_start=200, max_missed_end=50,
                      full_missed_start=50, full_missed_end=30):
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
            # opposite strand not allowed!
            if (cID == r.qID or (r.strand == '-' and same_strand_only)):
                yield HitItem(qID=r.qID, cID=cID)
                continue

            # regardless if whether is full-length (is_FL)
            # the query MUST be mapped fully (based on full_missed_start/end)
            if r.qStart > full_missed_start or (r.qLength-r.qEnd) > full_missed_end:
                yield HitItem(qID=r.qID, cID=cID)

            # full-length case: allow up to max_missed_start bp of 5' not aligned
            # and max_missed_end bp of 3' not aligned
            # non-full-length case: not really tested...don't use
            if is_FL and not alignment_missed_start_end_less_than_threshold(r,\
                            max_missed_start, max_missed_end, full_missed_start, full_missed_end):
                yield HitItem(qID=r.qID, cID=cID)
            else:
                cigar_str, ece_arr = eval_blasr_alignment(
                    record=r,
                    qver_get_func=qver_get_func,
                    qvmean_get_func=qvmean_get_func,
                    sID_starts_with_c=sID_starts_with_c,
                    qv_prob_threshold=qv_prob_threshold)

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

def daligner_against_ref2(query_dazz_handler, target_dazz_handler, la4ice_filename,
                         is_FL, sID_starts_with_c,
                         qver_get_func, qvmean_get_func, qv_prob_threshold=.03,
                         ece_penalty=1, ece_min_len=20, same_strand_only=True, no_qv_or_aln_checking=False,
                         max_missed_start=200, max_missed_end=50,
                         full_missed_start=50, full_missed_end=30):
    """
    Excluding criteria:
    (1) self hit
    (2) opposite strand hit  (should already be in the same orientation;
        can override with <same_strand_only> set to False)
    (3) less than 90% aligned or more than 50 bp missed

    Parameters:
      query_dazz_handler - query dazz handler in DalignRunner
      target_dazz_handler - target dazz handler in DalignRunner
      la4ice_filename - la4ice output of DalignRunner
      qver_get_func - returns a list of qvs of (read, qvname)
                      e.g. basQV.basQVcacher.get() or .get_smoothed()
      qvmean_get_func - which returns mean QV of (read, qvname)
    """
    for r in LA4IceReader(la4ice_filename):
        missed_q = r.qStart + r.qLength - r.qEnd
        missed_t = r.sStart + r.sLength - r.sEnd

        r.qID = query_dazz_handler[r.qID].split(' ')[0]
        r.sID = target_dazz_handler[r.sID].split(' ')[0]

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
        if (cID == r.qID or (r.strand == '-' and same_strand_only)):
            yield HitItem(qID=r.qID, cID=cID)
            continue

        # regardless if whether is full-length (is_FL)
        # the query MUST be mapped fully (based on full_missed_start/end)
        print "r.qStart:", r.qID, r.sID, r.qStart, full_missed_start, (r.qLength-r.qEnd), full_missed_end, r.qStart > full_missed_start or (r.qLength-r.qEnd) > full_missed_end

        if r.qStart > full_missed_start or (r.qLength-r.qEnd) > full_missed_end:
            yield HitItem(qID=r.qID, cID=cID)
            continue

        # this is used for partial_uc/nFL reads only
        # simply accepts hits from daligner for the nFL partial hits
        # testing shows that it does not affect much the Quiver consensus calling
        if no_qv_or_aln_checking:
            yield HitItem(qID=r.qID, cID=cID,
                          qStart=r.qStart, qEnd=r.qEnd,
                          missed_q=missed_q * 1. / r.qLength,
                          missed_t=missed_t * 1. / r.sLength,
                          fakecigar=1,
                          ece_arr=1)
            continue

        # full-length case: allow up to 200bp of 5' not aligned
        # and 50bp of 3' not aligned
        if (is_FL and not alignment_missed_start_end_less_than_threshold(r, \
             max_missed_start, max_missed_end, full_missed_start, full_missed_end)):
            yield HitItem(qID=r.qID, cID=cID)
        else:
            cigar_str, ece_arr = eval_blasr_alignment(
                record=r,
                qver_get_func=qver_get_func,
                sID_starts_with_c=sID_starts_with_c,
                qv_prob_threshold=qv_prob_threshold,
                qvmean_get_func=qvmean_get_func)
            #else: # don't use QV, just look at alignment

            if alignment_has_large_nonmatch(ece_arr, ece_penalty, ece_min_len):
                yield HitItem(qID=r.qID, cID=cID)
            else:
                yield HitItem(qID=r.qID, cID=cID,
                              qStart=r.qStart, qEnd=r.qEnd,
                              missed_q=missed_q * 1. / r.qLength,
                              missed_t=missed_t * 1. / r.sLength,
                              fakecigar=cigar_str,
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