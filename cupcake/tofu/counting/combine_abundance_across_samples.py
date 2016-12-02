__author__ = 'etseng@pacb.com'


import os, sys
from collections import defaultdict
from cupcake.tofu import compare_junctions
from cupcake.io import GFF
from bx.intervals import IntervalTree
from bx.intervals.cluster import ClusterTree

class MegaPBTree:
    """
    Structure for maintaining a non-redundant set of gene annotations
    Used to combine with different collapsed GFFs from different samples
    """
    def __init__(self, gff_filename, group_filename, internal_fuzzy_max_dist=0, self_prefix=None):
        self.gff_filename = gff_filename
        self.group_filename = group_filename
        self.self_prefix = self_prefix
        self.internal_fuzzy_max_dist = internal_fuzzy_max_dist
        self.record_d = dict((r.seqid, r) for r in GFF.collapseGFFReader(gff_filename))
        self.tree = defaultdict(lambda: {'+':IntervalTree(), '-':IntervalTree()}) # chr --> strand --> tree

        #print >> sys.stderr, "self.internal_fuzzy_max_dist is", internal_fuzzy_max_dist
        #raw_input()
        self.read_gff_as_interval_tree()
        self.group_info = MegaPBTree.read_group(self.group_filename, self.self_prefix) # ex: PB.1.1 --> [ RatHeart|i3_c123.... ]

    def read_gff_as_interval_tree(self):
        """
        Read a collapsed GFF file into an IntervalTree
        """
        for r in GFF.collapseGFFReader(self.gff_filename):
            self.tree[r.chr][r.strand].insert(r.start, r.end, r)

    @staticmethod
    def read_group(group_filename, group_prefix):
        group_info = {}
        with open(group_filename) as f:
            for line in f:
                pbid, members = line.strip().split('\t')
                if group_prefix is None:
                    group_info[pbid] = [x for x in members.split(',')]
                else:
                    group_info[pbid] = [group_prefix+'|'+x for x in members.split(',')]
        return group_info

    def match_record_to_tree(self, r):
        """
        r --- GMAPRecord
        tree --- dict of chromosome --> strand --> IntervalTree

        If exact match (every exon junction), return the matching GMAPRecord
        Otherwise return None
        *NOTE*: the tree should be non-redundant so can return as soon as exact match is found!
        """
        matches = self.tree[r.chr][r.strand].find(r.start, r.end)
        for r2 in matches:
            r.segments = r.ref_exons
            r2.segments = r2.ref_exons
            if compare_junctions.compare_junctions(r, r2, internal_fuzzy_max_dist=self.internal_fuzzy_max_dist) == 'exact': # is a match!
                return r2
        return None

    def add_sample(self, gff_filename, group_filename, sample_prefix, output_prefix):
        combined = [] # list of (r1 if r2 is None | r2 if r1 is None | longer of r1 or r2 if both not None)
        unmatched_recs = self.record_d.keys()

        for r in GFF.collapseGFFReader(gff_filename):
            match_rec = self.match_record_to_tree(r)
            if match_rec is not None:  # found a match! put longer of r1/r2 in
                combined.append((match_rec, r))
                try:
                    unmatched_recs.remove(match_rec.seqid)
                except ValueError:
                    pass # already deleted, OK, this happens for single-exon transcripts
            else:  # r is not present in current tree
                combined.append((None, r))
        # put whatever is left from the tree in
        for seqid in unmatched_recs:
            combined.append((self.record_d[seqid], None))

        # create a ClusterTree to re-calc the loci/transcripts
        final_tree = defaultdict(lambda: {'+': ClusterTree(0, 0), '-':ClusterTree(0, 0)})
        for i,(r1,r2) in enumerate(combined):
            if r2 is None or (r1 is not None and r1.end-r1.start > r2.end-r2.start):
                final_tree[r1.chr][r1.strand].insert(r1.start, r1.end, i)
            else:
                final_tree[r2.chr][r2.strand].insert(r2.start, r2.end, i)

        self.write_cluster_tree_as_gff(final_tree, combined, group_filename, sample_prefix, output_prefix)


    def write_cluster_tree_as_gff(self, cluster_tree, rec_list, group_filename2, sample_prefix2, output_prefix):
        """
        Write ClusterTree (chr --> dict --> (start, end, rec_list_index)) as collapsedGFF format
        Returns --- a new group_info!!!
        """
        group_info2 = MegaPBTree.read_group(group_filename2, sample_prefix2)
        new_group_info = {}
        f_out = open(output_prefix+'.gff', 'w')
        f_group = open(output_prefix+'.group.txt', 'w')
        f_mgroup = open(output_prefix + '.mega_info.txt', 'w')
        f_mgroup.write("pbid\t{0}\t{1}\n".format(self.self_prefix, sample_prefix2))
        loci_index = 0
        chroms = cluster_tree.keys()
        chroms.sort()
        for k in chroms:
            for strand in ('+', '-'):
                for _s, _e, rec_indices in cluster_tree[k][strand].getregions():
                    loci_index += 1
                    isoform_index = 0
                    for i in rec_indices:
                        isoform_index += 1
                        tID = "PB.{i}.{j}".format(i=loci_index,j=isoform_index)
                        r1, r2 = rec_list[i]
                        if r1 is None: # r2 is not None
                            r = r2
                            new_group_info[tID] = group_info2[r2.seqid]
                            f_mgroup.write("{tID}\tNA\t{group}\n".format(tID=tID, group=r2.seqid))
                        elif r2 is None: # r1 is not None
                            r = r1
                            new_group_info[tID] = self.group_info[r1.seqid]
                            f_mgroup.write("{tID}\t{group}\tNA\n".format(tID=tID, group=r1.seqid))
                        else: # both r1, r2 are not empty
                            r = r1 if (r1.end-r1.start > r2.end-r2.start) else r2
                            new_group_info[tID] = self.group_info[r1.seqid] + group_info2[r2.seqid]
                            f_mgroup.write("{tID}\t{group1}\t{group2}\n".format(tID=tID, group1=r1.seqid, group2=r2.seqid))

                        f_group.write("{tID}\t{members}\n".format(tID=tID, members=",".join(new_group_info[tID])))
                        f_out.write("{chr}\tPacBio\ttranscript\t{s}\t{e}\t.\t{strand}\t.\tgene_id \"PB.{i}\"; transcript_id \"{tID}\";\n".format(\
                            chr=k, s=r.start+1, e=r.end, strand=strand, tID=tID, i=loci_index))
                        for exon in r.ref_exons:
                            f_out.write("{chr}\tPacBio\texon\t{s}\t{e}\t.\t{strand}\t.\tgene_id \"PB.{i}\"; transcript_id \"{tID}\";\n".format(\
                                chr=k, s=exon.start+1, e=exon.end, strand=strand, tID=tID, i=loci_index))
        f_out.close()
        f_group.close()
        f_mgroup.close()
        return new_group_info


