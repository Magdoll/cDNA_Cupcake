__author__ = 'etseng@pacb.com'

import fire

from Bio import SeqIO
from cupcake.io.BLASRRecord import BLASRM5Reader
from cupcake.ice.ice_align_core import eval_blasr_alignment, alignment_has_large_nonmatch


#id      length  is_fl   stat    pbid

def is_a_hit(r, is_fl, max_missed_qstart=50, max_missed_qend=50, max_missed_5=400, max_missed_3=50, ece_penalty=1, ece_min_len=20):
    """
    r --- BLASRm5Record

    """
    if r.strand!='+': return False
    if r.qStart > max_missed_qstart: return False
    if r.qLength-r.qEnd > max_missed_qend: return False
    if is_fl:
        if r.sStart > max_missed_5: return False
        if r.sLength-r.sEnd > max_missed_3: return False
    cigar_str, ece = eval_blasr_alignment(r, lambda a,b,c: 1, lambda a,b: 1, False, 0)
    return not alignment_has_large_nonmatch(ece, ece_penalty, ece_min_len)


def process_blasr_hits(blasr_filename, is_fl, max_missed_qstart=50, max_missed_qend=50, max_missed_5=400, max_missed_3=50, ece_penalty=1, ece_min_len=20):
    """
    Note: blasr file must be grouped by qID!
    For FL reads -- whichever is the last fulfilling record gets picked
    For nFL reads --- record all instances
    """
    last_qID = None
    records_seen = []
    for r in BLASRM5Reader(blasr_filename):
        if r.qID != last_qID: # write out current record
            if last_qID is not None and len(records_seen) > 0:
                yield records_seen
            last_qID = r.qID
            records_seen = []
        if is_fl:
            if is_a_hit(r, True, max_missed_qstart, max_missed_qend, max_missed_5, max_missed_3, ece_penalty, ece_min_len):
                records_seen = [r]
        else: # nFL
            if is_a_hit(r, False, max_missed_qstart, max_missed_qend, max_missed_5, max_missed_3, ece_penalty, ece_min_len):
                records_seen.append(r)

def process_blasr_file(fasta_filename, blasr_filename, is_fl, output_filename, is_pbid=True, output_mode='w', max_missed_qstart=50, max_missed_qend=50, max_missed_5=400, max_missed_3=50, ece_penalty=1, ece_min_len=20):
    unaligned = dict((r.id,len(r.seq)) for r in SeqIO.parse(open(fasta_filename), 'fasta'))

    with open(output_filename, output_mode) as f:
        f.write("id\tlength\tis_fl\tstat\tpbid\n")
        for records in process_blasr_hits(blasr_filename, is_fl, max_missed_qstart, max_missed_qend, max_missed_5, max_missed_3, ece_penalty, ece_min_len):
            stat = 'unique' if len(records) == 1 else "ambiguous"
            for r in records:
                if is_pbid:
                    pbid = r.sID.split('|')[0]
                else:
                    pbid = r.sID
                f.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(\
                    r.qID, r.qLength, 'Y' if is_fl else 'N', stat, pbid))
            del unaligned[records[0].qID]

        for _id,_len in unaligned.iteritems():
            f.write("{0}\t{1}\t{2}\tNA\tNA\n".format(_id, _len, 'Y' if is_fl else 'N'))

if __name__ == "__main__":
    fire.Fire(process_blasr_file, name='process')

