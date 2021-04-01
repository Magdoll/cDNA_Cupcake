
"""
Call variant based on a list of SAMMPileUpRecord where list[i] is the record of i-th position

Most of the code follows Juliet's code at
https://github.com/PacificBiosciences/minorseq/blob/develop/src/AminoAcidCaller.cpp

"""
import os, sys
import scipy.stats as stats
from collections import Counter, namedtuple

BHtuple = namedtuple('BHtuple', ['pval', 'record'])

class MPileUPVariant(object):
    def __init__(self, record_list, min_cov, err_sub, expected_strand, pval_cutoff=0.01, bhFDR=None):
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
        self.bhFDR       = bhFDR   # is None, this is not used; other wise do Benjamini–Hochberg
        self.expected_strand = expected_strand


        self.prep_records()
        self.positions_to_call = self.get_positions_to_call()

        # must first call positions to call, then prep records, then number of tests
        self.number_of_tests = sum(self.record_by_pos[pos].clean_type for pos in self.positions_to_call)

        self.variant = {} # position --> in sorted order, (base, count)
        self.ref_base = {} # position --> ref base

        self.call_variant()


    def is_in_or_near_HP(self, pos, hp_size=4):
        """
        We define a HP region as stretches of 4 or more same nucleotides
        :return: True/False for in/hear HP region
        """
        def find_hp_region_size(cur):
            if cur not in self.record_by_pos: return 0
            end = cur+1
            while end in self.record_by_pos and self.record_by_pos[end].ref == self.record_by_pos[cur].ref:
                end += 1
            start = cur-1
            while start in self.record_by_pos and self.record_by_pos[start].ref == self.record_by_pos[cur].ref:
                start -= 1
            # hp region is from start+1 to end
            return end-(start+1)

        return (find_hp_region_size(pos) >= hp_size) or \
               (find_hp_region_size(pos-1) >= hp_size) or \
               (find_hp_region_size(pos+1) >= hp_size)


    def get_positions_to_call(self):
        """
        Identify list of positions to try to call SNPs. Must have:
        1. minimum coverage >= min_cov
        2. the first and second most frequent base are NOT an indel
        3. not next to or inside a homopolymer region
        4. has at least two or more keys
        """
        positions_to_call = []
        for pos in self.record_by_pos:
            if self.record_by_pos[pos].clean_type < 2: continue # only one base at this position, skip
            elif self.record_by_pos[pos].clean_cov < self.min_cov: continue # insufficient cov, skip
            else:
                # find the first and second most freq base in the "non-clean" counts
                m = self.record_by_pos[pos].clean_counts.most_common()
                # ex: m = [('a', 10), ('-ct', 20), ('+t', 10)]
                if m[0][0][0]in ('+','-') or m[1][0][0] in ('+','-') or self.is_in_or_near_HP(pos): continue
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
            if self.expected_strand == '+-':
                # for metagenomics, we don't care the strand
                # so instead we will convert everything to upper case later in the counts
                bases = 'ATCG'
            elif self.expected_strand == '+':
                bases = 'ATCG'
            elif self.expected_strand == '-':
                bases = 'atcg'

            if self.expected_strand == '+-':
                # convert lower case to upper case
                for k in 'atcg':
                    if k in r.counts:
                        r.counts[k.upper()] += r.counts[k]
                        del r.counts[k]

            r.clean_counts = Counter(r.counts)
            keys = list(r.counts.keys())
            for k in keys:
                if k not in bases:
                    del r.clean_counts[k]
            r.clean_cov = sum(r.clean_counts.values())
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
        if self.bhFDR is None: # use Bonferroni correction
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
        else: # use Benjamini–Hochberg procedure
            # see: https://www.statisticshowto.com/benjamini-hochberg-procedure/
            pval_dict = {} # (pos, base) -> BHtuple(pval, record)
            for pos in self.positions_to_call:
                r = self.record_by_pos[pos]
                for base, count in r.clean_counts.most_common()[1:]:
                    assert not base.startswith('+') and not base.startswith('-') # clean counts should NOT have indels
                    exp = r.clean_cov * self.err_sub
                    odds, pval = stats.fisher_exact([[count, r.clean_cov-count], [exp, r.clean_cov-exp]], alternative='greater')
                    if pval <= self.pval_cutoff:      # With this filtration, the sequencing errors position will not be stored in pval_dict.
                        pval_dict[(pos, base)] = BHtuple(pval=pval, record=r)

            # now we have all the pvals, rank them
            keys_pos_base = list(pval_dict.keys())
            keys_pos_base.sort(key=lambda x: pval_dict[x].pval)
            self.number_of_tests = len(keys_pos_base)
            # find the largest p value that is smaller than the critical value.
            largest_good_rank1 = 0
            for rank0,(pos, base) in enumerate(keys_pos_base):
                pval = pval_dict[(pos, base)].pval
                bh_val = ((rank0+1)/self.number_of_tests) * self.bhFDR # Only significant positions will be used to adjust bh_val
                if pval < bh_val:
                    largest_good_rank1 = rank0+1
                    print(f"pos:{pos} base:{base} pval:{pval} bh:{bh_val}")
            for (pos,base) in keys_pos_base[:largest_good_rank1]:
                r = pval_dict[(pos,base)].record
                if pos not in self.variant:
                    self.ref_base[pos] = r.ref
                    self.variant[pos] = [r.clean_counts.most_common()[0]]
                self.variant[pos] += [(base, r.clean_counts[base])]


class MagMPileUPVariant(MPileUPVariant):
    def __init__(self, record_list, min_cov, err_sub, expected_strand, pval_cutoff=0.01, bhFDR=None):
        self.ref_name = {} # position --> ref contig
        super().__init__(record_list, min_cov, err_sub, expected_strand, pval_cutoff, bhFDR)

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
        if self.bhFDR is None: # use Bonferroni correction
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
                    self.ref_name[pos] = r.chr

        else: # use Benjamini–Hochberg procedure
            # see: https://www.statisticshowto.com/benjamini-hochberg-procedure/
            pval_dict = {} # (pos, base) -> BHtuple(pval, record)
            for pos in self.positions_to_call:
                r = self.record_by_pos[pos]
                for base, count in r.clean_counts.most_common()[1:]:
                    assert not base.startswith('+') and not base.startswith('-') # clean counts should NOT have indels
                    exp = r.clean_cov * self.err_sub
                    odds, pval = stats.fisher_exact([[count, r.clean_cov-count], [exp, r.clean_cov-exp]], alternative='greater')
                    if pval <= self.pval_cutoff:      # With this filtration, the sequencing errors position will not be stored in pval_dict.
                        pval_dict[(pos, base)] = BHtuple(pval=pval, record=r)
            # now we have all the pvals, rank them
            keys_pos_base = list(pval_dict.keys())
            keys_pos_base.sort(key=lambda x: pval_dict[x].pval)
            self.number_of_tests = len(keys_pos_base)
            # find the largest p value that is smaller than the critical value.
            largest_good_rank1 = 0
            for rank0,(pos, base) in enumerate(keys_pos_base):
                pval = pval_dict[(pos, base)].pval
                bh_val = ((rank0+1)/self.number_of_tests) * self.bhFDR # Only significant positions will be used to adjust bh_val
                if pval < bh_val:
                    largest_good_rank1 = rank0+1
                    print(f"pos:{pos} base:{base} pval:{pval} bh:{bh_val}")
            for (pos,base) in keys_pos_base[:largest_good_rank1]:
                r = pval_dict[(pos,base)].record
                if pos not in self.variant:
                    self.ref_base[pos] = r.ref
                    self.ref_name[pos] = r.chr
                    self.variant[pos] = [r.clean_counts.most_common()[0]]
                self.variant[pos] += [(base, r.clean_counts[base])]




