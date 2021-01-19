__author__ = 'etseng@pacb.com'

"""
ex:
?CAACTTATTGTCCC 0   1   0   0
?CACCT?ATTGTCCC 0   1   0   0
?CA?CTTATTGTCCC 0   1   0   0
G?ACTGCGGAACTTT 0   1   0   0
GTAAT?CG?A?TTCT 0   1   0   0
GTAATGCG?AACTTT 0   2   0   0

we have:
1. list of haplotypes (ex: "AAT", "A?T", "ATG") and the weight of each haplotype
2. list of SNPs at each position
3. user input of ploidity N

start with assumption of 2 allele:
how do we find the center
"""
import pdb
import os, sys
import itertools
import numpy as np

import networkx as nx
from collections import Counter
from phasing.io.VariantPhaser import Haplotypes


def make_haplotype_counts(isoform_tally):
    """
    :param hap_obj: Haplotype object
    :param isoform_tally: list of (isoform, dict of haplotype count), ex: {'PB.45.1': {0:10, 1:20}}
    :return: Counter of haplotype index --> total count (of FL reads)
    """
    hap_count = Counter() # haplotype index --> total count
    for tally in isoform_tally.values():
        for hap_index, count in tally.items():
            hap_count[hap_index] += count
    return hap_count


def calc_hap_diff(haplotype_strings, cur_hap):
    """
    :param haplotype_strings: list of haplotype strings
    :param cur_hap: hap string to compare against

    :return: np.array of diffs against <cur_hap>
    """
    return np.array([sum(cur_hap[i]!=s for i,s in enumerate(hap_str)) for hap_str in haplotype_strings])


def get_hap_model(haplotype_strings):
    hap_len = len(haplotype_strings[0])
    tally = [Counter() for pos in range(hap_len-1)]
    for s in haplotype_strings:
        for pos in range(hap_len-1):
            if s[pos]!='?' and s[pos+1]!='?':
                tally[pos][(s[pos], s[pos+1])] += 1
    return tally


def infer_haplotypes_via_exhaustive_diploid_only(hap_obj, variants):
    """
    ONLY WORKS FOR DIPLOID RIGHT NOW!!
    This is an alternative to the more heuristic infer_haplotypes_via_min_diff()
    Brute force search through all solutions which is 2^N (N=number of variants)

    :param haplotype_obj: haplotype object to be error corrected
    :param variants: length N list of variants which we will use to permutate

    1. permute through all possible diploid pairs
    2. calc diff of each haplotype string against the diploid hap strings
    3. find the (hap1,hap2) pair that minimizes sum of diffs

    :return: sum_of_diff array, hap_count_ordered (the counts are fake since not used later)
    """
    # sanity check that every position has exactly two possibilities
    #for x in variants:
    #    if len(x)!=2: raise Exception, "variants must be diploid!"

    haplotype_strings = hap_obj.haplotypes
    nonpartial_haps = [s for s in haplotype_strings if all(x!='?' for x in s)]
    nonpartial_haps_already_tried = set()

    if len(nonpartial_haps) == 0:
        return None, None

    #is_list = False
    #hap_len = len(haplotype_strings[0])
    #
    #if len(nonpartial_haps) == 0: # special case! there are NO non-partial haplotype strings!
    #    freq_model = get_hap_model(haplotype_strings)
    #    print freq_model
    #    nonpartial_haps = itertools.product(*variants)
    #    is_list = True

    best_diff_arr, best_sum_of_diff = None, np.inf
    best_hap1, best_hap2 = None, None
    for hap1 in nonpartial_haps:
        if hap1 in nonpartial_haps_already_tried: continue
        #if is_list:
        #    is_likely = True
        #    # see how likely hap1 is given the freq model
        #    for pos in xrange(hap_len-1):
        #        if (hap1[pos], hap1[pos+1]) not in freq_model[pos]:
        #            print (hap1[pos], hap1[pos + 1]), "not in pos", pos
        #            is_likely = False
        #            break
        #    if is_likely:
        #        hap1 = "".join(hap1)
        #        print >> sys.stderr, "Trying hap1: ", hap1
        #    else:
        #        continue
        # hap2 is the other choice since it's diploid
        hap2 = ''
        for i in range(len(variants)):
            if hap1[i]==variants[i][0]: hap2 += variants[i][1]
            else: hap2 += variants[i][0]
        # now calculate the sum diffs
        diff_arr = np.array([calc_hap_diff(haplotype_strings, hap1),
                             calc_hap_diff(haplotype_strings, hap2)])
        sum_of_diff = diff_arr.min(axis=0).sum()
        if sum_of_diff < best_sum_of_diff:
            best_diff_arr, best_sum_of_diff = diff_arr, sum_of_diff
            best_hap1, best_hap2 = hap1, hap2
        nonpartial_haps_already_tried.add(hap1)
        nonpartial_haps_already_tried.add(hap2)

    best_hap1_index, msg1 = hap_obj.match_or_add_haplotype(best_hap1)
    best_hap2_index, msg2 = hap_obj.match_or_add_haplotype(best_hap2)
    return best_diff_arr, [(best_hap1_index,-1), (best_hap2_index,-2)]


def infer_haplotypes_via_min_diff(haplotype_strings, hap_count, ploidy, max_diff, min_percent):
    """
    :param haplotype_strings: list of haplotype strings
    :param hap_count: Counter object of hap_index --> count
    :param ploidy: user given ploidy (max number of true haplotypes)

    1. sort haplotypes based on counts. select hap0 to be most frequent one.
    2. calc diff of each haplotype against hap0. calc sum of diffs.
    3. add in hap1 as second most frequent, calc sum of diffs.
    4. ...repeat until sum of diffs gets worse, or ploidy reached, or the next haplotype is less than <min_percent> reads .
    """
    #pdb.set_trace()
    # we are to ignore all partial haps
    partial_haps = [i for i in range(len(haplotype_strings)) if any(s=='?' for s in haplotype_strings[i])]
    for k in partial_haps:
        del hap_count[k]

    if len(hap_count) == 0:
        return None, None

    hap_count_ordered = hap_count.most_common() # now sorted in desc (hap_index, count)
    hap_count_total = sum([count for hap_index,count in hap_count_ordered])
    cur_hap = haplotype_strings[hap_count_ordered[0][0]]
    len_hap = len(cur_hap)

    if len_hap < 2*max_diff:
        max_diff = 1

    diff_arr = np.array([calc_hap_diff(haplotype_strings, cur_hap)])
    sum_of_diff = diff_arr.sum()
    #print sum_of_diff
    for cur_hap_i, cur_count in hap_count_ordered[1:ploidy]:
        new_diff_arr = np.append(diff_arr, [calc_hap_diff(haplotype_strings, haplotype_strings[cur_hap_i])], axis=0)
        new_sum_of_diff = new_diff_arr.min(axis=0).sum()
        #print new_sum_of_diff
        if cur_count < hap_count_total * min_percent:
            break
        if (sum_of_diff - new_sum_of_diff) < max_diff:
            break
        else:
            sum_of_diff = new_sum_of_diff
            diff_arr = new_diff_arr

    return diff_arr, hap_count_ordered


def error_correct_haplotypes(hap_obj, isoform_tally, diff_arr, hap_count_ordered):

    # create new hap_obj and old_to_new_map dict
    new_hap_obj = Haplotypes(hap_obj.hap_var_positions, hap_obj.ref_at_pos, hap_obj.count_of_vars_by_pos)
    old_to_new_map = {}
    for i, j in enumerate(diff_arr.argmin(axis=0)):
        # haplotype i maps to haplotype hap_count_ordered[j][0]
        k = hap_count_ordered[j][0]
        new_hap_index, msg = new_hap_obj.match_or_add_haplotype(hap_obj.haplotypes[k])
        old_to_new_map[i] = new_hap_index

    # now create a new isoform_tally
    new_isoform_tally = {}
    for k,v in isoform_tally.items():
        new_isoform_tally[k] = Counter()
        for old_hap_index, count in v.items():
            if old_hap_index not in old_to_new_map:
                print("Discarding: {0}".format(hap_obj.haplotypes[old_hap_index]), file=sys.stderr)
                continue
            new_hap_index = old_to_new_map[old_hap_index]
            new_isoform_tally[k][new_hap_index] += count
    return old_to_new_map, new_hap_obj, new_isoform_tally


# Liz: BELOW ARE OBSOLETE.
# def make_haplotype_counts(isoform_tally):
#     """
#     :param hap_obj: Haplotype object
#     :param isoform_tally: list of (isoform, dict of haplotype count), ex: {'PB.45.1': {0:10, 1:20}}
#     :return: dict of haplotype index --> total count (of FL reads)
#     """
#     hap_count = Counter() # haplotype index --> total count
#     for tally in isoform_tally.itervalues():
#         for hap_index, count in tally.iteritems():
#             hap_count[hap_index] += count
#     return hap_count
#
#
# def make_haplotype_counts_into_graph(hap_obj, isoform_tally):
#     """
#     :param hap_obj: Haplotype object
#     :param isoform_tally: list of (isoform, dict of haplotype count), ex: {'PB.45.1': {0:10, 1:20}}
#     """
#     G = nx.DiGraph()
#     for tally in isoform_tally.itervalues():
#         for hap_index, count in tally.iteritems():
#             hap_str = hap_obj.haplotypes[hap_index]
#             n = len(hap_str)
#             for i in xrange(n-1):
#                 # skip any non-variant-assigned locations
#                 if hap_str[i]=='?' or hap_str[i+1]=='?': continue
#                 node1 = (i, hap_str[i])
#                 node2 = (i+1, hap_str[i+1])
#                 if node1 in G and node2 in G[node1]:
#                     G[node1][node2]['weight'] += count
#                 else:
#                     G.add_edge(node1, node2, weight=count)
#     return G
#
# def enumerate_allele_candidates(G, cur_node, cur_str, other_str):
#     if G.out_degree(cur_node) == 0: # reached a sink
#         print cur_str
#     else:
#         es = G[cur_node].items()
#         es.sort(key=lambda x: x[1]['weight'], reverse=True)
#         # ex: [((2, 'A'), {'weight': 23}), ((2, 'C'), {'weight': 1})]
#         for (node, _dict) in es:
#             enumerate_allele_candidates(G, node, cur_str+node[1])
#
#
# def make_haplotype_graph_nonpartial_only(haplotype_strings, err_sub, max_diff_allowed):
#     """
#     :param haplotype_strings: list of haplotype strings (in order), ex: ["AATT", "TCGG", "A?TT"]
#     :param err_sub: given substitution error rate
#     :param max_diff_allowed: maximium base diff allowed to connect the two strings (graph edges)
#     :return: a nx.Graph where nodes are haplotype indices and connecting edges indicate string similarity
#
#     This is a prep function for creating a similarity graph where we can use to eventually identify the true alleles.
#     *NOTE:* Hap strings with missing bases ('?') are ignored.
#     """
#     G = nx.Graph()
#     n = len(haplotype_strings)
#     hap_str_len = len(haplotype_strings[0])
#     if hap_str_len > 2*max_diff_allowed:
#         max_diff = max(max_diff_allowed, err_sub * hap_str_len)
#     else:
#         max_diff = err_sub * hap_str_len
#     partial_haps = filter(lambda i: any(s=='?' for s in haplotype_strings[i]), xrange(n))
#     for i1 in xrange(n):
#         if i1 in partial_haps: continue # skip ones with '?'
#         G.add_node(i1)
#         for i2 in xrange(i1+1, n):
#             if i2 in partial_haps: continue # skip ones with '?'
#             s1 = haplotype_strings[i1]
#             s2 = haplotype_strings[i2]
#             sim = sum((s1[k]==s2[k]) for k in xrange(hap_str_len))
#             if (hap_str_len-sim) < max_diff:
#                 G.add_edge(i1, i2, weight=sim)
#     return G, partial_haps
#
# def error_correct_haplotypes(G, partial_haps, hap_count, hap_obj, isoform_tally):
#     """
#     :param G: Graph of haplotype str similarities (via make_haplotype_graph_nonpartial_only), using only nonpartial hap strings
#     :param partial_haps: list of indices of haps that are partial (have '?') and not in G. need to be imputed.
#     :param hap_count: dict of haplotype index --> total count (via make_haplotype_counts)
#     :param hap_obj: Haplotype object
#     :param isoform_tally: list of (isoform, dict of haplotype count), ex: {'PB.45.1': {0:10, 1:20}}
#
#     :return: new_to_map_mapping, new_hap_obj, new_isoform_tally
#
#     1. first identify cliques from G
#     2. each clique becomes a new haplotype; pick the most common one are the representative
#     3. error correct the other clique members to become this rep (adding counts to it)
#     """
#     cliques = [comm for comm in nx.k_clique_communities(G, 2)]
#     # adding all orphans as a clique by itself
#     nodes_left = set(G.nodes())
#
#     new_hap_obj = Haplotypes(hap_obj.hap_var_positions, hap_obj.ref_at_pos, hap_obj.count_of_vars_by_pos)
#     # for each clique, pick the one that has the most counts
#     old_to_new_map = {}
#     for i,members in enumerate(cliques):
#         for hap_index in members: nodes_left.remove(hap_index)
#         stuff = [(hap_index, hap_count[hap_index]) for hap_index in members]
#         stuff.sort(key=lambda x: x[1], reverse=True)
#         # choose the most abundant haplotype (none of them have '?')
#         chosen_hap_index = stuff[0][0] # init to first one
#         # --- below not needed anymore becuz we're only using nonpartial strings with no '?' ---
#         #for hap_index, count_ignore in stuff:
#         #    hap_str = hap_obj.haplotypes[hap_index]
#         #    if all(s!='?' for s in hap_str):
#         #        chosen_hap_index = hap_index
#         #        break
#         new_hap_index, msg = new_hap_obj.match_or_add_haplotype(hap_obj.haplotypes[chosen_hap_index])
#         for x in members: old_to_new_map[x] = new_hap_index
#
#     # look through all leftover nodes, all of which do not have '?'
#     for hap_index in nodes_left:
#         new_hap_index, msg = new_hap_obj.match_or_add_haplotype(hap_obj.haplotypes[hap_index])
#         old_to_new_map[hap_index] = new_hap_index
#
#     # impute partial hap strings, assign them if they are sufficiently close to a hap
#     for hap_index in partial_haps:
#         # min_score = 3 means the partial hap must have at least two base matches with the highest scoring one
#         new_hap_index = new_hap_obj.impute_haplotype(hap_obj.haplotypes[hap_index], min_score=3)
#         if new_hap_index is not None:
#             old_to_new_map[hap_index] = new_hap_index
#         else: # could not impute, discard, so do nothing
#             pass
#
#     # now create a new isoform_tally
#     new_isoform_tally = {}
#     for k,v in isoform_tally.iteritems():
#         new_isoform_tally[k] = Counter()
#         for old_hap_index, count in v.iteritems():
#             if old_hap_index not in old_to_new_map:
#                 print >> sys.stderr, "Discarding: {0}".format(hap_obj.haplotypes[old_hap_index])
#                 continue
#             new_hap_index = old_to_new_map[old_hap_index]
#             new_isoform_tally[k][new_hap_index] += count
#     return old_to_new_map, new_hap_obj, new_isoform_tally



