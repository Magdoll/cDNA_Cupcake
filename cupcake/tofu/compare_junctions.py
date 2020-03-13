#!/usr/bin/env python
__author__ = 'etseng@pacb.com'

import pdb

def overlaps(s1, s2):
    return max(0, min(s1.end, s2.end) - max(s1.start, s2.start))

def compare_junctions(r1, r2, internal_fuzzy_max_dist=0, max_5_diff=999999, max_3_diff=999999):
    """
    r1, r2 should both be BioReaders.GMAPSAMRecord

    super
    exact
    subset
    partial
    nomatch

    <internal_fuzzy_max_dist> allows for very small amounts of diff between internal exons
    useful for chimeric & slightly bad mappings
    """
    #assert internal_fuzzy_max_dist >= 0
    #assert max_5_diff >= 0
    #assert max_3_diff >= 0
    strand = r1.strand
    found_overlap = False
    # super/partial --- i > 0, j = 0
    # exact/partial --- i = 0, j = 0
    # subset/partial --- i = 0, j > 0
    for i,x in enumerate(r1.segments):
        # find the first matching r2, which could be further downstream
        for j,y in enumerate(r2.segments):
            if i > 0 and j > 0: break
            if overlaps(x, y) > 0:
                found_overlap = True
                break
        if found_overlap:
            break
    if not found_overlap: return "nomatch"
    # now we have r1[i] matched to r2[j]
    # check that r1[i] and r2[j] match within 5'/3' max diff
    if strand == '+':
        if abs(r1.segments[i].start-r2.segments[j].start) > max_5_diff or \
           abs(r1.segments[i].end-r2.segments[j].end) > max_3_diff:
           return "partial"
    else:
        if abs(r1.segments[i].start-r2.segments[j].start) > max_3_diff or \
           abs(r1.segments[i].end-r2.segments[j].end) > max_5_diff:
            return "partial"

    # if just one exon, then must have less than 5'/3' diff
    #pdb.set_trace()


    if len(r1.segments) == 1:
        if len(r2.segments) == 1:
            # if both are single-exon, check that they have less than 5'/3' diff
            if abs(r1.segments[0].start-r2.segments[0].start) <= max_5_diff and \
               abs(r1.segments[0].end-r2.segments[0].end) <= max_3_diff:
                return "exact"
            else:
                return "partial"
        else:
            # case: r1 single exon, r2 multi-exon
            # check that the matching exon is within difference
            # and that the r1 exon does NOT overlap with previous or later exons in r2
            if abs(r1.segments[0].end-r2.segments[j].end) <= max_3_diff and \
                abs(r1.segments[0].start-r2.segments[j].start) <= max_5_diff and \
                    ((j==0 and r1.segments[0].end < r2.segments[1].start) or
                      (j>=1 and r1.segments[0].start > r2.segments[j-1].end and (j==len(r2.segments)-1 or r1.segments[0].end < r2.segments[j+1].start))):
                return "subset"
            else:
                return "partial"
    else:
        if len(r2.segments) == 1:
            # case: r1 multi exon, r2 single exon
            # r1.segments[i] matches r2.segments[0]
            # need to check that the r2, being single exon, did not overlap with previous or later r1 exons
            # and also the matching r1, r2 exon does not have huge 5'/3' difference
            if (i==0 or r1.segments[i-1].end < r2.segments[0].start) and \
               (i==len(r1.segments)-1 or r1.segments[i].start > r2.segments[0].end) and \
               (abs(r1.segments[i].start-r2.segments[0].start) <= max_5_diff) and \
               (abs(r1.segments[i].end-r2.segments[0].end) <= max_3_diff):
                return "super"
            else:
                return "partial"
        else: # both r1 and r2 are multi-exon, check that all remaining junctions agree
            k = 0
            while i+k+1 < len(r1.segments) and j+k+1 < len(r2.segments):
                if abs(r1.segments[i+k].end-r2.segments[j+k].end)>internal_fuzzy_max_dist or \
                   abs(r1.segments[i+k+1].start-r2.segments[j+k+1].start)>internal_fuzzy_max_dist:
                    return "partial"
                k += 1
            #print i, j, k
            # check that the last matching exon, the ends are with max 5'/3' diff
            if strand == '+':
                if abs(r1.segments[i+k].end-r2.segments[j+k].end) > max_3_diff:
                    return "partial"
            else:
                if abs(r1.segments[i+k].end-r2.segments[j+k].end) > max_5_diff:
                    return "partial"

            if i+k+1 == len(r1.segments):
                if j+k+1 == len(r2.segments):
                    if i == 0:
                        if j == 0: return "exact"
                        else: return "subset"    # j > 0
                    else: return "super"
                else: # r1 is at end, r2 not at end
                    if i == 0: return "subset"
                    else:  # i > 0
                        if abs(r1.segments[i+k-1].end-r2.segments[j+k-1].end)>internal_fuzzy_max_dist or \
                           abs(r1.segments[i+k].start-r2.segments[j+k].start)>internal_fuzzy_max_dist:
                            return "partial"
                        else:
                            return "concordant"
            else: # r1 not at end, r2 must be at end
                if j == 0: return "super"
                else:
                    if abs(r1.segments[i+k-1].end-r2.segments[j+k-1].end)>internal_fuzzy_max_dist or \
                        abs(r1.segments[i+k].start-r2.segments[j+k].start)>internal_fuzzy_max_dist:
                        return "partial"
                    else:
                        return "concordant"


