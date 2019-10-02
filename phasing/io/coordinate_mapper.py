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


def make_exons_from_base_mapping(mapping, start, end, strand):
    """
    mapping is 0-based index on transcript  --> 0-based index  on genome
    however beware of strand!
    """

    output = [mapping[start]]
    for i in range(start+1, end):
        cur_pos, cur_is_junction= mapping[i]
        if cur_is_junction and mapping[i]!=output[-1]:
            # if the last position is the same, DON'T APPEND (was an indel)
            output.append(mapping[i])
    cur_pos, cur_is_junction = mapping[end]
    if mapping[end]!=output[-1]:
        output.append(mapping[end])

    # remember for Interval it is 0-based start, 1-based end
    # if len(output) is odd, must be 1bp into an exon
    # ex: [(xxx,True), (xxx,True), (xxx,False)] or
    #     [.....(xxx,True), xxx(True)]
    #print output
    if len(output)==1:
        output = [output[0], output[0]] # just duplicate it
    elif len(output)%2==1:
        if output[0][1] and output[1][1]:
            output.insert(0, output[0])
        elif output[-1][1] and output[-2][1]:
            output.append(output[-1])
    #    print "modified:", output
    if strand == '+':
        return [Interval(output[i][0],output[i+1][0]+1) for i in range(0, len(output), 2)]
    else: # - strand
        return [Interval(output[i][0],output[i-1][0]+1)  for i in range(len(output)-1,-1,-2)]



def get_base_to_base_mapping_from_sam(exons, cigar_string, qStart, qEnd, strand, include_junction_info=False):
    """
    For PacBio data which can have indels w.r.t genome =___=

    ex:
        cigar: 1S105M407N548M
        sStart-sEnd: 948851-949911
        qStart-qEnd: 2-655
        segments: [Interval(start=948851, end=948956), Interval(start=949363, end=949911)]

    Returns: dict of 0-based position --> 0-based ref position
    """
    cur_exon_i = 0
    cur_nt_loc = qStart
    cur_genome_loc = exons[0].start

    start_soft_clip = qStart > 0

    last_base_is_junction = False
    qLen = qEnd

    mapping = {}

    for num, s in iter_cigar_string(cigar_string):
        if s == 'S': # soft clipping at the ends, ignore
            if start_soft_clip:
                assert num == qStart
                for i in range(num): mapping[i] = (cur_genome_loc, False)
                start_soft_clip = False
            else:
                # soft clipping at the end
                # advance the mapping but not cur_nt_loc (otherwise will be diff from qEnd)
                for i in range(num): 
                    mapping[cur_nt_loc+i] = (cur_genome_loc, False)
                    #cur_nt_loc += 1
                    #print cur_nt_loc
                    # for soft clipping at the end, do NOT progress cur_nt_loc!
                    # we are now "outside" the alignment, otherwise
                    # assert cur_nt_loc == qEnd will be wrong at the end
                #cur_nt_loc -= 1
                qLen += num # query length must be qEnd + soft clipped end
        elif s == 'N': # intron, move to next ref exon
            mapping[cur_nt_loc-1] = (mapping[cur_nt_loc-1][0], True)
            assert cur_genome_loc == exons[cur_exon_i].end
            cur_exon_i += 1
            cur_genome_loc = exons[cur_exon_i].start
            last_base_is_junction = True
        elif s == 'M':
            # for the next "num" matches are all 1:1
            for i in range(num):
                if cur_nt_loc in mapping and mapping[cur_nt_loc][1]:
                    # if this is true, then last mapping must be 'D' and was a junction
                    # so we do nothing -- keep it
                    pass
                else:
                    mapping[cur_nt_loc] = (cur_genome_loc, last_base_is_junction)
                last_base_is_junction = False
                cur_nt_loc += 1
                cur_genome_loc += 1
            assert cur_genome_loc <= exons[cur_exon_i].end
        elif s == 'I': # insertion w.r.t to genome
            for i in range(num):
                mapping[cur_nt_loc] = (cur_genome_loc, last_base_is_junction)
                cur_nt_loc += 1
                last_base_is_junction = False
        elif s == 'D': # deletion w.r.t. to genome
            # if last_base_is_junction is True, we want to make sure it makes it in mapping
            mapping[cur_nt_loc] = (cur_genome_loc, last_base_is_junction)
            last_base_is_junction = False
            cur_genome_loc += num
# BELOW IS WRONG
#            for i in xrange(num):
#                mapping[cur_nt_loc] = cur_genome_loc
#                cur_genome_loc += 1
            assert cur_genome_loc <= exons[cur_exon_i].end
    assert cur_nt_loc == qEnd or (cur_nt_loc==qEnd-1 and s=='S')

    if strand == '-':
        mapping = dict((qLen-1-k, v) for k,v in mapping.items())

    if not include_junction_info:
        for k in mapping:
            mapping[k] = mapping[k][0]

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

