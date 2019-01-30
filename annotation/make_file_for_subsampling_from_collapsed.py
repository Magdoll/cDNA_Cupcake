#!/usr/bin/env python
import os, sys
from csv import DictReader
from Bio import SeqIO


def make_file_for_subsample(input_prefix, output_filename, matchAnnot_parsed=None, sqanti_class=None, include_single_exons=False):
    """
    Two files must exist: .abundance.txt and .rep.fq so we can make the length
    """
    count_filename = input_prefix + '.abundance.txt'
    fq_filename = input_prefix + '.rep.fq'

    if not include_single_exons:
        from cupcake.io.GFF import collapseGFFReader
        gff_filename = input_prefix + '.gff'
        print >> sys.stderr, "Reading {0} to exclude single exons...".format(gff_filename)
        good_ids = []
        for r in collapseGFFReader(gff_filename):
            if len(r.ref_exons) >= 2:
                good_ids.append(r.seqid)

    if not os.path.exists(count_filename):
        print >> sys.stderr, "Cannot find {0}. Abort!".format(count_filename)
        sys.exit(-1)

    if not os.path.exists(fq_filename):
        print >> sys.stderr, "Cannot find {0}. Abort!".format(fq_filename)
        sys.exit(-1)

    if matchAnnot_parsed is not None and not os.path.exists(matchAnnot_parsed):
        print >> sys.stderr, "Cannot find {0}. Abort!".format(matchAnnot_parsed)
        sys.exit(-1)

    if sqanti_class is not None and not os.path.exists(sqanti_class):
        print >> sys.stderr, "Cannot find {0}. Abort!".format(sqanti_class)
        sys.exit(-1)

    if matchAnnot_parsed is not None:
        match_dict = dict((r['pbid'],r) for r in DictReader(open(matchAnnot_parsed), delimiter='\t'))
    elif sqanti_class is not None:
        print >> sys.stderr, "Reading {0} to get gene/isoform assignment...".format(sqanti_class)
        match_dict = {}
        for r in DictReader(open(sqanti_class), delimiter='\t'):
            if r['associated_transcript'] == 'novel':
                refisoform = 'novel_'+r['isoform']
            else:
                refisoform = r['associated_transcript']
            match_dict[r['isoform']] = {'refgene': r['associated_gene'],
                                     'refisoform': refisoform}
    else:
        match_dict = None

    seqlen_dict = dict((r.id.split('|')[0],len(r.seq)) for r in SeqIO.parse(open(fq_filename),'fastq'))
    
    h = open(output_filename, 'w')
    if matchAnnot_parsed is None and sqanti_class is None:
        h.write("pbid\tpbgene\tlength\tfl_count\n")
    else:
        h.write("pbid\tpbgene\tlength\trefisoform\trefgene\tfl_count\n")
    f = open(count_filename)
    while True:
        cur = f.tell()
        if not f.readline().startswith('#'):
            f.seek(cur)
            break
    for r in DictReader(f, delimiter='\t'):
        if not include_single_exons and r['pbid'] not in good_ids:
            print >> sys.stderr, "Exclude {0} because single exon.".format(r['pbid'])
            continue
        h.write("{0}\t{1}\t{2}\t".format(r['pbid'], r['pbid'].split('.')[1], seqlen_dict[r['pbid']]))
        if matchAnnot_parsed is not None or sqanti_class is not None:
            m = match_dict[r['pbid']]
            h.write("{0}\t{1}\t".format(m['refisoform'], m['refgene']))
        h.write("{0}\n".format(r['count_fl']))
    h.close()

    print >> sys.stderr, "Output written to {0}.".format(output_filename)

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser("Make subsample-ready file from ToFU/Iso-Seq collapsed output")
    parser.add_argument("-i", "--input_prefix", default="hq_isoforms.fastq.no5merge.collapsed.min_fl_2.filtered", help="Collapsed prefix (default: hq_isoforms.fastq.no5merge.collapsed.min_fl_2.filtered)")
    parser.add_argument("-o", "--output_filename", default="hq_isoforms.fastq.no5merge.collapsed.min_fl_2.filtered.for_subsampling.txt", help="Output filename (default: hq_isoforms.fastq.no5merge.collapsed.min_fl_2.filtered.for_subsampling.txt")
    parser.add_argument("-m1", "--matchAnnot_parsed", default=None, help="MatchAnnot parsed output (default: None)")
    parser.add_argument("-m2", "--sqanti_class", default=None, help="SQANTI classification file (default: None)")
    parser.add_argument("--include_single_exons", default=False, action="store_true", help="Include single exons (default: OFF)")

    args = parser.parse_args()
    make_file_for_subsample(args.input_prefix, args.output_filename, args.matchAnnot_parsed, args.sqanti_class, args.include_single_exons)



