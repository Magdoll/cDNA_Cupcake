#!/usr/bin/env python
import os, sys
from csv import DictReader
from Bio import SeqIO

def make_file_for_subsample(input_prefix, output_filename, matchAnnot_parsed=None):
    """
    Two files must exist: .abundance.txt and .rep.fq so we can make the length
    """
    count_filename = input_prefix + '.abundance.txt'
    fq_filename = input_prefix + '.rep.fq'

    if not os.path.exists(count_filename):
        print >> sys.stderr, "Cannot find {0}. Abort!".format(count_filename)
        sys.exit(-1)

    if not os.path.exists(fq_filename):
        print >> sys.stderr, "Cannot find {0}. Abort!".format(fq_filename)
        sys.exit(-1)

    if matchAnnot_parsed is not None and not os.path.exists(matchAnnot_parsed):
        print >> sys.stderr, "Cannot find {0}. Abort!".format(matchAnnot_parsed)
        sys.exit(-1)

    if matchAnnot_parsed is None:
        match_dict = None
    else:
        match_dict = dict((r['pbid'],r) for r in DictReader(open(matchAnnot_parsed), delimiter='\t'))

    seqlen_dict = dict((r.id.split('|')[0],len(r.seq)) for r in SeqIO.parse(open(fq_filename),'fastq'))
    
    h = open(output_filename, 'w')
    if matchAnnot_parsed is None:
        h.write("pbid\tpbgene\tlength\tfl_count\n")
    else:
        h.write("pbid\tpbgene\tlength\trefisoform\trefgene\tfl_count\n")
    f = open(count_filename)
    for i in xrange(14): assert f.readline().startswith('#') # read the header section
    for r in DictReader(f, delimiter='\t'):
        h.write("{0}\t{1}\t{2}\t".format(r['pbid'], r['pbid'].split('.')[1], seqlen_dict[r['pbid']]))
        if matchAnnot_parsed is not None:
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
    parser.add_argument("-m", "--matchAnnot_parsed", default=None, help="MatchAnnot parsed output (default: None)")

    args = parser.parse_args()
    make_file_for_subsample(args.input_prefix, args.output_filename, args.matchAnnot_parsed)



