
"""
Call variant based on a list of SAMMPileUpRecord where list[i] is the record of i-th position

Most of the code follows Juliet's code at
https://github.com/PacificBiosciences/minorseq/blob/develop/src/AminoAcidCaller.cpp

"""
import os, sys
import scipy.stats as stats
from collections import Counter

class MPileUPVariant(object):
    def __init__(self, record_list, min_cov, err_sub, expected_strand, pval_cutoff=0.01):
        """
        :param record_list: list of SAMMPileUpRecord
        :param min_cov: minimum coverage to call variant
        :param err_sub: substitution error, right now a fixed float
        :param expected_strand: expected strand of the transcript (+ or -)
        """
        self.record_by_pos = dict((r.pos, r) for r in record_list)
        self.min_cov = min_cov
        self.err_sub = err_sub
        self.pval_cutoff = pval_cutoff
        self.expected_strand = expected_strand

        self.prep_records()
        self.positions_to_call = self.get_positions_to_call()

        # must first call positions to call, then prep records, then number of tests
        self.number_of_tests = sum(self.record_by_pos[pos].clean_type for pos in self.positions_to_call)

        self.variant = {} # position --> in sorted order, (base, count)
        self.ref_base = {} # position --> ref base

        self.call_variant()


    def get_positions_to_call(self):
        """
        Identify list of positions to try to call SNPs. Must have:
        1. minimum coverage >= min_cov
        2. the first or second most frequent base is NOT an indel
        3. has at least two or more keys
        """
        positions_to_call = []
        for pos in self.record_by_pos:
            if self.record_by_pos[pos].clean_type < 2: continue # only one base at this position, skip
            elif self.record_by_pos[pos].clean_cov < self.min_cov: continue # insufficient cov, skip
            else:
                # find the first and second most freq base in the "non-clean" counts
                m = self.record_by_pos[pos].clean_counts.most_common()
                # ex: m = [('a', 10), ('-ct', 20), ('+t', 10)]
                if m[0][0][0]in ('+','-') or m[1][0][0] in ('+','-'): continue
                else:
                    positions_to_call.append(pos)
        return positions_to_call

    def prep_records(self):
        """
        Prepare the records by:
        1. remove all 'N' bases
        2. remove all bases that were not on the expected strand
        3. remove all indels

        Creates three new vars: clean_counts, clean_cov, clean_type
        DOES NOT ALTER the original counts or other variables!!!

        If + strand, then ATCG
        If - strand, then atcg
        """
        for pos in self.record_by_pos:
            r = self.record_by_pos[pos]
            if self.expected_strand == '+':
                bases = 'ATCG'
            else:
                bases = 'atcg'

            r.clean_counts = Counter(r.counts)
            keys = r.counts.keys()
            for k in keys:
                if k not in bases:
                    del r.clean_counts[k]
            r.clean_cov = sum(r.clean_counts.itervalues())
            r.clean_type = len(r.clean_counts)

    def call_variant(self):
        """
        mirrors AminoAcidCaller::CallVariants() in
        https://github.com/PacificBiosciences/minorseq/blob/develop/src/AminoAcidCaller.cpp

        For each position (that has sufficient coverage),
         do Fisher exact test w/ correction
         if p-val < threshold, then store it.

        Stores results in self.variant as:

        self.variant[position] = desc list of (base, count).
        NOTE: base must be either all in lower case (which means - strand)
              or call upper case (+ strand).
              If - strand and ('a', 10), it means the ref base in A on the + strand,
              and the transcript should be T on the - strand.

        Only positions with more than the ref base is stored.
        """
        for pos in self.positions_to_call:
            r = self.record_by_pos[pos]
            alt_variant = []
            for base, count in r.clean_counts.most_common()[1:]:
                assert not base.startswith('+') and not base.startswith('-') # clean counts should NOT have indels
                exp = r.clean_cov * self.err_sub
                odds, pval = stats.fisher_exact([[count, r.clean_cov-count], [exp, r.clean_cov-exp]], alternative='greater')
                pval *= self.number_of_tests
                if pval < self.pval_cutoff: # store variant if below cutoff
                    alt_variant.append((base, count))
            if len(alt_variant) > 0: # only record this variant if there's at least two haps
                self.variant[pos] = [r.clean_counts.most_common()[0]] + alt_variant
                self.ref_base[pos] = r.ref







