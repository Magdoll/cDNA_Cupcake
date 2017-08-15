#!/usr/bin/env python
__author__ = 'etseng@pacb.com'

import os, sys
from bx.intervals import IntervalTree, Interval
from Bio import SeqIO

try:
    import BED
except ImportError:
    print >> sys.stderr, "Cannot find BED.py! Please make sure you have cDNA_Cupcake/sequence in $PYTHONPATH."
try:
    import BioReaders
except ImportError:
    print >> sys.stderr, "Cannot find BioReaders.py! Please make sure you have cDNA_Cupcake/sequence in $PYTHONPATH."


def get_probe_hit(tree, gene_info, r):
    """
    Given a dict tree (from read_probe_bed) and a GMAP SAM record
    Go through each exon and find probes that hit it

    Return: (number of probes hit), (total number of bases overlapping with probes), (genes seen)
    """
    probes_seen = set()
    genes_seen = set()
    base_hit = 0
    if r.sID not in tree: return 0, 0, set()
    for e in r.segments:
        hits = tree[r.sID].find(e.start, e.end)
        if len(hits) == 0: continue
        for i,intl in hits:
            probes_seen.add(i)
            genes_seen.add(gene_info[i])
            base_hit += min(e.end,intl.end)-max(e.start,intl.start)
    return len(probes_seen), base_hit, genes_seen

def read_probe_bed(bed_filename, start_base=0, end_base=1):
    """
    Read a probe BED file <chrom>, <start>, <end>
    Return dict of chrom --> IntervalTree w/ data=(index, interval)
    """
    tree = {}
    gene_info = {}
    i = 0
    reader = BED.SimpleBEDReader(bed_filename, start_base, end_base)
    for r in reader:
        if r.chr not in tree: tree[r.chr] = IntervalTree()
        tree[r.chr].add(r.start, r.end, (i,Interval(r.start, r.end)))
        if r.name is not None:
            gene_info[i] = r.name
        i += 1
    return tree, gene_info

def calc_ontarget_rate(tree, gene_info, input_fasta, sam_filename, output_filename=None):

    query_len_dict = dict((r.id, len(r.seq)) for r in SeqIO.parse(open(input_fasta),'fasta'))

    if output_filename is None:
        f = sys.stdout
    else:
        f = open(output_filename, 'w')
    f.write("read_id\tread_len\tis_fl\tnum_probe\tnum_base_overlap\tloci\tgenes\n")
    reader = BioReaders.GMAPSAMReader(sam_filename, True,query_len_dict=query_len_dict)
    for r in reader:
        if r.sID=='*': continue
        num_probe,base_hit,genes_seen=get_probe_hit(tree, gene_info, r)
        f.write("{0}\t{1}\tY\t{2}\t{3}\t{4}:{5}-{6}\t{7}\n".format(r.qID,r.qLen,num_probe,base_hit,r.sID,r.sStart,r.sEnd,",".join(genes_seen)))

    f.close()

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser("Calculate Probe Hit from SAM alignment + probe BED")
    parser.add_argument("bed_filename")
    parser.add_argument("input_fasta")
    parser.add_argument("sam_filename")
    parser.add_argument("--start_base", choices=['0', '1'], required=True, help="Start base is 0 or 1-based index")
    parser.add_argument("--end_base", choices=['0', '1'], required=True, help="End base is 0 or 1-based index")
    parser.add_argument("-o", "--output", help="Output filename (default: stdout)")

    args = parser.parse_args()

    tree, gene_info = read_probe_bed(args.bed_filename, int(args.start_base), int(args.end_base))

    calc_ontarget_rate(tree, gene_info, args.input_fasta, args.sam_filename, args.output)
