__author__ = 'etseng@pacb.com'

import os, sys
import re
from bx.intervals import Interval

def sanity_check(seq, motifs, tally):
    for motif in motifs:
        n = len(motif)
        for pos in tally[motif]:
            assert seq[pos:pos+n] == motif


def count_repeats_for_motif(seq, motif, tally, intervals=None):
    """
    seq --- plain sequence to search for the repeats (motifs)
    motif --- plain sequence of repeat, ex: CGG, AGG
    intervals --- 0-based start, 1-based end of Intervals to search motif in
    """
    if intervals is None: # use the whole sequence
        intervals = [Interval(0, len(seq))]

    new_intl = []
    for intl in intervals:
        cur = seq[intl.start:intl.end]
        prev_end = intl.start
        found_flag = False
        for m in re.finditer(motif, cur):
            tally[motif].append(intl.start + m.start())
            if m.start() > prev_end:
                # new interval is prev_end (0-based), m.start() (1-based)
                new_intl.append(Interval(prev_end, intl.start + m.start()))
            prev_end = intl.start + m.end()
            found_flag = True
        if not found_flag:
            new_intl.append(intl)
    return new_intl


def count_repeats(seq, motifs):

    tally = {}
    intervals = None
    for motif in motifs:
        tally[motif] = []
        intervals = count_repeats_for_motif(seq, motif, tally, intervals)

    return intervals, tally


def mask_repeats_in_sequence(seq, intervals, mask='_'):
    if len(intervals) == 0: return mask * len(seq)
    newseq = ''
    newseq += mask * intervals[0].start

    for i in xrange(len(intervals)-1):
        newseq += seq[intervals[i].start:intervals[i].end]
        newseq += mask * (intervals[i+1].start-intervals[i].end)

    newseq += seq[intervals[-1].start:intervals[-1].end]
    newseq += mask * (len(seq)-intervals[-1].end)
    return "".join(newseq)


