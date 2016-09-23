#!/usr/bin/env python

"""match_w_annotation.py: functions for categorizing read alignments to annotations. mainly for bacteria."""

__author__ = "Elizabeth Tseng"
__copyright__ = "Copyright 2016, cDNA_Cupcake"
__email__ = "etseng@pacb.com"

import os, sys
import pdb
from collections import defaultdict, namedtuple
from csv import DictReader
from Bio import SeqIO
from bx.intervals import IntervalTree
from bx.intervals.cluster import ClusterTree
import BioReaders  # should already have cDNA_Cupcake/sequence in $PYTHONPATH to get this

AMatch = namedtuple('AMatch', 'name strand start end record')

def calc_overlap(s1, e1, s2, e2):
    return min(e1, e2) - max(s1, s2)

def calc_overlap_ratio1(s1, e1, s2, e2):
    return (min(e1, e2)-max(s1, s2))*1./(e1-s1)

def check_multigene(overlaps, min_overlap_bp=200, min_query_overlap=.5, min_gene_overlap=.8):
    """
    overlaps is a list of: (gene, overlap_bp, overlap_gene_ratio, overlap_query_ratio)
    """
    if all(x[1]>=min_overlap_bp or x[2]>=min_gene_overlap for x in overlaps):
        new_name = "poly-"+overlaps[0][0]+'-'+overlaps[1][0]
        return new_name
    elif overlaps[0][3] >= min_query_overlap: # first gene covers 50% of query
        return overlaps[0][0]
    elif overlaps[1][3] >= min_query_overlap:
        return overlaps[1][0]
    else:
        return "novel"

def match_w_annotation(t, r, info, min_overlap_bp=200, min_query_overlap=.5, min_gene_overlap=.8):
    """
    Input:
       t -- dict of chr -> strand -> bx.intervals.IntervalTree
       r -- GMAP SAM record
       info --- dict of (start, end, gene_name)

    Return: AMatch(matching_name, strand, start, end, record)
    matching_name --- "novel" if it matches 0 genes; 
        gene_name if it matches that 1 gene;
        gene1_gene2_... if it matches multiple genes
    """
    s, e = r.sStart, r.sEnd
    matches = t[r.sID][r.flag.strand].find(s, e)
    if len(matches) == 0: return AMatch("novel", r.flag.strand, s, e, r)
    elif len(matches) == 1: 
        s2, e2, gene2 = info[matches[0]]
        # check that >= 50% of the query overlaps
        if calc_overlap_ratio1(s, e, s2, e2) >= min_query_overlap:
            return AMatch(gene2, r.flag.strand, s, e, r)
        else:
            return AMatch("novel", r.flag.strand, s, e, r)
    else: # matches 2+ genes
        #pdb.set_trace()
        overlaps = [] # list of (gene, overlap_bp, overlap_gene_ratio, overlap_query_ratio)
        for gene2 in matches:
            s2, e2, gene2 = info[gene2]
            o = calc_overlap(s, e, s2, e2)
            overlaps.append((gene2, o, o*1./(e2-s2), o*1./(e-s)))
        
        if len(overlaps) == 2: # must cover both genes by 200
            return AMatch(check_multigene(overlaps, min_overlap_bp, min_query_overlap, min_gene_overlap), r.flag.strand, s, e, r)
        else: # 3+ genes
            flag = check_multigene(overlaps, min_overlap_bp, min_query_overlap, min_gene_overlap)
            if flag!="novel": return AMatch(flag, r.flag.strand, s, e, r)
            # start shaving off genes at ends to see if they now match
            for i in xrange(1, len(overlaps)-1):
                flag = check_multigene(overlaps[i:], min_overlap_bp, min_query_overlap, min_gene_overlap)
                if flag!="novel": return AMatch(flag, r.flag.strand, s, e, r)
            for i in xrange(len(overlaps)-1, 1, -1):
                flag = check_multigene(overlaps[:i], min_overlap_bp, min_query_overlap, min_gene_overlap)
                if flag!="novel": return AMatch(flag, r.flag.strand, s, e, r)
            return AMatch("novel", r.flag.strand, s, e, r)


def categorize_aln_by_annotation(gene_annotation_file, input_fasta, input_sam, output_prefix, min_overlap_bp=200, min_query_overlap=.5, min_gene_overlap=.8):

    t = defaultdict(lambda: {'+':IntervalTree(), '-':IntervalTree()})  # chr -> strand -> IntervalTree
    info = {}

    #reader = DictReader(open('ProteinTable149_154224.txt'),delimiter='\t')
    for r in DictReader(open(gene_annotation_file),delimiter='\t'):
        if r['#Replicon Name']!='chr':
            print >> sys.stderr, "Ignore", r
            continue
        info[r['Locus tag']] = (int(r['Start']), int(r['Stop']), r['Locus tag'])
        t[r['Replicon Accession']][r['Strand']].add(int(r['Start']), int(r['Stop']), r['Locus tag'])

    #pdb.set_trace()

    result = defaultdict(lambda: []) # gene -> list of rec
    d = dict((r.id,len(r.seq)) for r in SeqIO.parse(open(input_fasta),'fasta'))

    reader = BioReaders.GMAPSAMReader(input_sam, True, query_len_dict=d)
    for r in reader:
        ans = match_w_annotation(t, r, info, min_overlap_bp, min_query_overlap, min_gene_overlap)
        # ans is AMatch(name, strand, start, end, record)
        result[ans.name].append(ans)

    novel_ct = defaultdict(lambda: {'+':ClusterTree(0,0), '-':ClusterTree(0,0)})
    novel_list = []
    novel_index = 0

    f = open(output_prefix + '.sam', 'w')
    f.write(reader.header)
    f1 = open(output_prefix + '.report.txt', 'w')
    f1.write("id\tread_group\tgene_name\tserial_number\tstrand\tstart\tend\n")
    for k,v in result.iteritems():
        # v is: list of AMatch(name, strand, start, end, record)
        if k=='novel':
            # write novel later, we are grouping them by loci first
            #tagRG='novel'
            for x in v:
                novel_ct[x.record.sID][x.strand].insert(x.start, x.end, novel_index)
                novel_index += 1
                novel_list.append(x)
            continue
        elif k.startswith('poly-'): tagRG='poly'
        else: tagRG='single'
        v.sort(key=lambda x: (x.start, x.end), reverse=True if v[0].strand=='-' else False) # sort by start, then end
        for i,x in enumerate(v):
            f.write("{0}\tSN:Z:{1:06d}\tRG:Z:{2}\tgn:Z:{3}\n".format(x.record.record_line, i+1, tagRG, k))
            if x.strand == '+':
                f1.write("{0}\t{1}\t{2}\t{3:06d}\t{4}\t{5}\t{6}\n".format(\
                    x.record.qID, tagRG, k, i+1, x.strand, x.start+1, x.end))
            else: # - strand, start is end, end is start
                f1.write("{0}\t{1}\t{2}\t{3:06d}\t{4}\t{5}\t{6}\n".format(\
                    x.record.qID, tagRG, k, i+1, x.strand, x.end, x.start+1))

    # now write the novel stuff, grouped by regions
    novel_region_index = 1
    for d1 in novel_ct.itervalues():
        for ct in d1.itervalues():
            gn = 'novel-' + str(novel_region_index)
            for _start, _end, _indices in ct.getregions():
                v = [novel_list[ind] for ind in _indices]
                v.sort(key=lambda x: (x.start, x.end), reverse=True if v[0].strand=='-' else False) # sort by start, then end
                for i,x in enumerate(v):
                    f.write("{0}\tSN:Z:{1:06d}\tRG:Z:{2}\tgn:Z:{3}\n".format(x.record.record_line, i+1, "novel", gn))
                    if x.strand == '+':
                        f1.write("{0}\t{1}\t{2}\t{3:06d}\t{4}\t{5}\t{6}\n".format(\
                            x.record.qID, "novel", gn, i+1, x.strand, x.start+1, x.end))
                    else:
                        f1.write("{0}\t{1}\t{2}\t{3:06d}\t{4}\t{5}\t{6}\n".format(\
                            x.record.qID, "novel", gn, i+1, x.strand, x.end, x.start+1))
                novel_region_index += 1

    f.close()
    f1.close()

    print >> sys.stderr, "Output written to:", f.name
    print >> sys.stderr, "Output written to:", f1.name


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser("Match alignment with annotation. Categorize and Report.")
    parser.add_argument("gene_annotation_file", help="Gene Annotation Text File")
    parser.add_argument("input_fasta", help="Input Fasta")
    parser.add_argument("input_sam", help="Input SAM")
    parser.add_argument("output_prefix", help="Output Prefix")
    parser.add_argument("--min_query_overlap", type=float, default=.5, help="Minimum query overlap (default: 0.5)")
    parser.add_argument("--min_poly_overlap_bp", type=int, default=200, help="Minimum gene overlap, in bp (default: 200 bp)")
    parser.add_argument("--min_poly_gene_overlap", type=float, default=.8, help="Minimum polycistronic gene overlap, in ratio (default: 0.8)")

    args = parser.parse_args()
    categorize_aln_by_annotation(args.gene_annotation_file, args.input_fasta, args.input_sam, args.output_prefix, \
                                 args.min_poly_overlap_bp, args.min_query_overlap, args.min_poly_gene_overlap)