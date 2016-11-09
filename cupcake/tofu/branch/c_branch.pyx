import numpy as np
cimport numpy as np
from cpython cimport bool
from cupcake.tofu.branch.intersection_unique import IntervalTreeUnique, Interval

DTYPE = np.int
ctypedef np.int_t DTYPE_t

def exon_finding(np.ndarray[DTYPE_t, ndim=1] baseC, np.ndarray[DTYPE_t, ndim=1] altC_neg, np.ndarray[DTYPE_t, ndim=1] altC_pos, np.ndarray[DTYPE_t, ndim=1] transC, int size, int threshSplit, int threshBase, int offset):
    cdef int i
    cdef bool tag = False
    cdef index = 0
    exon_tree = IntervalTreeUnique()
    for i in xrange(size):
        if (baseC[i] > threshBase and not tag):  # a new exon!
            #print "new exon found at position", i
            e_start = i
            tag = True
        if tag:
            if i == size - 1: # reached the end of genome, end of exon too
                exon_tree.insert_interval(Interval(e_start+offset, size+offset, index))
                index += 1
                tag = False
            elif baseC[i] <= threshBase: # end of exon at i-1
                exon_tree.insert_interval(Interval(e_start+offset, i+offset, index))
                index += 1
                tag = False
            elif baseC[i] > 0 and (altC_pos[i] > threshSplit or altC_neg[i+1] < -threshSplit): # alt. junction found!
                # end the current exon at i and start a new one at i + 1
                print "alt. junction found at", i
                exon_tree.insert_interval(Interval(e_start+offset, i+1+offset, index))
                index += 1
                e_start = i + 1
    return exon_tree

def intervals_all_adjacent(x):
    """
    x --- a list of Intervals
    Return True if all are adjacent

    ex:
    [Interval(2083850, 2083974),
     Interval(2083974, 2083983),
     Interval(2083983, 2083988),
     Interval(2083988, 2084011)]
    """
    cdef int i
    for i in xrange(len(x)-1):
        if x[i].end != x[i+1].start:
            return False
    return True

def exon_matching(exon_tree, ref_exon, match_extend_tolerate_left, match_extend_tolerate_right, intervals_adjacent=True):
    """
    exon_tree --- an IntervalTree made from .baseC/.altC using exon detection; probably only short read data
    ref_exon --- an Interval representing an exon; probably from PacBio
    match_extend_tolerate --- maximum difference between the matched start/end

    find a continuous exon path (consisting of 1 or more nodes for which the intervals must be adjacent)
    in exon_tree that matches to ref_exon
    """
    matches = exon_tree.find(ref_exon.start, ref_exon.end)
    if len(matches) == 0: # likely due to very low coverage on transcript
        return None
    # check that all matches are adjacent (no splicing! this just one integral exon)
    if (not intervals_adjacent) or intervals_all_adjacent(matches):
        # check if the ends differ a little, if so, extend to min/max
        for i in xrange(len(matches)):
            d_start = abs(matches[i].start - ref_exon.start)
            #print "matching {0} to {1}".format(matches[i].start, ref_exon.start)
            #pdb.set_trace()
            if d_start <= match_extend_tolerate_left: # now find the furthest end that satisfies the results
                for j in xrange(len(matches)-1, i-1, -1):
                    if abs(matches[j].end - ref_exon.end) <= match_extend_tolerate_right:
                        return matches[i:(j+1)]
        return None
    else: # ack! could not find evidence for this :<
        return None
              



                                       
        
