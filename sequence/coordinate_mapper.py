import os, sys
import pdb
import bisect
from bx.intervals import Interval
from Bio.Seq import Seq

def iter_cigar_string(cigar_string):
    num = cigar_string[0]
    for s in cigar_string[1:]:
        if str.isalpha(s):
            yield int(num), s
            num = ''
        else:
            num += s

def get_base_to_base_mapping_from_sam(exons, cigar_string, qStart, qEnd):
    """
    For PacBio data which can have indels w.r.t genome =___=

    ex:
        cigar: 1S105M407N548M
        sStart-sEnd: 948851-949911
        qStart-qEnd: 2-655
        segments: [Interval(start=948851, end=948956), Interval(start=949363, end=949911)]
    """
    cur_exon_i = 0
    cur_nt_loc = qStart
    cur_genome_loc = exons[0].start

    start_soft_clip = qStart > 0

    mapping = {}

    for num, s in iter_cigar_string(cigar_string):
        if s == 'S': # soft clipping at the ends, ignore
            if start_soft_clip:
                assert num == qStart
                for i in xrange(num): mapping[i] = cur_genome_loc
                start_soft_clip = False
            else:
                # soft clipping at the end
                # advance the mapping but not cur_nt_loc (otherwise will be diff from qEnd)
                for i in xrange(num): 
                    mapping[cur_nt_loc+i] = cur_genome_loc
                    #cur_nt_loc += 1
                    #print cur_nt_loc
                cur_nt_loc -= 1
        elif s == 'N': # intron, move to next ref exon
            assert cur_genome_loc == exons[cur_exon_i].end
            cur_exon_i += 1
            cur_genome_loc = exons[cur_exon_i].start
        elif s == 'M':
            # for the next "num" matches are all 1:1
            for i in xrange(num):
                mapping[cur_nt_loc] = cur_genome_loc
                cur_nt_loc += 1
                cur_genome_loc += 1
            assert cur_genome_loc <= exons[cur_exon_i].end
        elif s == 'I': # insertion w.r.t to genome
            for i in xrange(num):
                mapping[cur_nt_loc] = cur_genome_loc
                cur_nt_loc += 1
        elif s == 'D': # deletion w.r.t. to genome
            for i in xrange(num):
                mapping[cur_nt_loc] = cur_genome_loc
                cur_genome_loc += 1
            assert cur_genome_loc <= exons[cur_exon_i].end
    assert cur_nt_loc == qEnd

    return mapping


def get_exon_coordinates(exons, start, end):
    """
    Return the set of "exons" (genome location) that 
    is where the nucleotide start-end is

    start is 0-based
    end is 1-based
    exons is a set of Interval (0-based start, 1-based end)
    """
    acc_lens = [0] # ex: [0, 945, 1065, 1141, 1237] accumulative length of exons
    len_of_transcript = 0
    for e in exons:
        _len = e.end - e.start
        acc_lens.append(acc_lens[-1] + _len)
        len_of_transcript += _len
    # confirm that start-end is in the range of the transcript!
    assert 0 <= start < end <= len_of_transcript + 30 # allow a 30-bp slack due to PacBio indels

    end = min(end, len_of_transcript) # trim it to the end if necessary (for PacBio)


    i = bisect.bisect_right(acc_lens, start) 
    j = bisect.bisect_right(acc_lens, end) 

    # starts at i-th exon and ends at j-th exon, i and j are both 1-based
    # for the first exon, the offset is start-acc+e.start
    # for the last exon, the end point is end-acc+e.start
    if i == j:
        return [Interval(start-acc_lens[i-1]+exons[i-1].start, 
                end-acc_lens[i-1]+exons[i-1].start)]
    else:
        if j >= len(exons):  # the end is the end
            return [Interval(start-acc_lens[i-1]+exons[i-1].start, exons[i-1].end)] + \
                    exons[i:] 
        else:
            return [Interval(start-acc_lens[i-1]+exons[i-1].start, exons[i-1].end)] + \
                exons[i:j-1] + \
               [Interval(exons[j-1].start, end-acc_lens[j-1]+exons[j-1].start)]
    
def consistute_genome_seq_from_exons(genome_dict, _chr, exons, strand):
    """
    genome_dict is expected to be SeqReaders.LazyFastaReader
    exons is a list of [Interval(start, end)]
    """
    seq = ''
    genome_seq = genome_dict[_chr].seq
    for e in exons:
        seq += str(genome_seq[e.start:e.end])

    seq = Seq(seq)
    if strand == '+':
        return seq.tostring()
    else:
        return seq.reverse_complement().tostring()

