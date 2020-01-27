#!/usr/bin/env python
import os, sys
from csv import DictReader
from Bio import SeqIO

def read_demux_fl_count_file(filename):
    d = {}
    reader = DictReader(open(filename), delimiter=',')
    assert 'id' in reader.fieldnames
    for r in reader:
        d[r['id']] = r
    samples = reader.fieldnames
    samples.remove('id')
    return d, samples

def make_file_for_subsample(input_prefix, output_prefix, demux_file=None, matchAnnot_parsed=None, sqanti_class=None, include_single_exons=False):
    """
    Two files must exist: .abundance.txt and .rep.fq so we can make the length
    """
    count_filename = input_prefix + '.abundance.txt'

    rep_filenames = [(input_prefix + '.rep.fq', 'fastq'), (input_prefix + '.rep.fastq','fastq'), \
                     (input_prefix + '.rep.fa', 'fasta'), (input_prefix + '.rep.fasta','fasta')]

    rep_filename = None
    rep_type = None
    for x,type in rep_filenames:
        if os.path.exists(x):
            rep_filename = x
            rep_type = type

    if rep_filename is None:
        print("Expected to find input fasta or fastq files {0}.rep.fa or {0}.rep.fq. Not found. Abort!".format(input_prefix), file=sys.stderr)
        sys.exit(-1)

    if not include_single_exons:
        from cupcake.io.GFF import collapseGFFReader
        gff_filename = input_prefix + '.gff'
        print("Reading {0} to exclude single exons...".format(gff_filename), file=sys.stderr)
        good_ids = []
        for r in collapseGFFReader(gff_filename):
            if len(r.ref_exons) >= 2:
                good_ids.append(r.seqid)

    if demux_file is None and not os.path.exists(count_filename):
        print("Cannot find {0}. Abort!".format(count_filename), file=sys.stderr)
        sys.exit(-1)

    if matchAnnot_parsed is not None and not os.path.exists(matchAnnot_parsed):
        print("Cannot find {0}. Abort!".format(matchAnnot_parsed), file=sys.stderr)
        sys.exit(-1)

    if sqanti_class is not None and not os.path.exists(sqanti_class):
        print("Cannot find {0}. Abort!".format(sqanti_class), file=sys.stderr)
        sys.exit(-1)

    if matchAnnot_parsed is not None:
        match_dict = dict((r['pbid'],r) for r in DictReader(open(matchAnnot_parsed), delimiter='\t'))
        for k in match_dict: match_dict[k]['category'] = match_dict[k]['score']
    elif sqanti_class is not None:
        print("Reading {0} to get gene/isoform assignment...".format(sqanti_class), file=sys.stderr)
        match_dict = {}
        for r in DictReader(open(sqanti_class), delimiter='\t'):
            if r['associated_transcript'] == 'novel':
                refisoform = 'novel_'+r['isoform']
            else:
                refisoform = r['associated_transcript']
            match_dict[r['isoform']] = {'refgene': r['associated_gene'],
                                        'refisoform': refisoform,
                                        'category': r['structural_category']}
    else:
        match_dict = None

    seqlen_dict = dict((r.id.split('|')[0],len(r.seq)) for r in SeqIO.parse(open(rep_filename),rep_type))

    to_write = {}
    if demux_file is None:
        to_write['all'] = {}
        f = open(count_filename)
        while True:
            cur = f.tell()
            if not f.readline().startswith('#'):
                f.seek(cur)
                break
        for r in DictReader(f, delimiter='\t'):
            if r['pbid'] in good_ids or include_single_exons:
                to_write['all'][r['pbid']] = r['count_fl']
    else:
        d, samples = read_demux_fl_count_file(demux_file)
        for s in samples: to_write[s] = {}
        for pbid, d2 in d.items():
            for s in samples:
                if pbid in good_ids or include_single_exons:
                    to_write[s][pbid] = d2[s]

    for sample in to_write:
        h = open(output_prefix+'.'+sample+'.txt', 'w')
        if matchAnnot_parsed is None and sqanti_class is None:
            h.write("pbid\tpbgene\tlength\tfl_count\n")
        else:
            h.write("pbid\tpbgene\tlength\trefisoform\trefgene\tcategory\tfl_count\n")
        for pbid in to_write[sample]:
            if matchAnnot_parsed is not None or sqanti_class is not None:
                if pbid not in match_dict:
                    print("Ignoring {0} because not on annotation (SQANTI/MatchAnnot) file.".format(pbid), file=sys.stdout)
                    continue
                m = match_dict[pbid]
                h.write("{0}\t{1}\t{2}\t".format(pbid, pbid.split('.')[1], seqlen_dict[pbid]))
                h.write("{0}\t{1}\t{2}\t".format(m['refisoform'], m['refgene'], m['category']))
            else:
                h.write("{0}\t{1}\t{2}\t".format(pbid, pbid.split('.')[1], seqlen_dict[pbid]))
            h.write("{0}\n".format(to_write[sample][pbid]))
        h.close()
        print("Output written to {0}.".format(h.name), file=sys.stderr)

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser("Make subsample-ready file from Iso-Seq collapsed output")
    parser.add_argument("-i", "--input_prefix", default="hq_isoforms.fastq.no5merge.collapsed.min_fl_2.filtered", help="Collapsed prefix (default: hq_isoforms.fastq.no5merge.collapsed.min_fl_2.filtered)")
    parser.add_argument("-o", "--output_prefix", default="output.for_subsampling", help="Output prefix (default: output.for_subsampling")
    parser.add_argument("-m1", "--matchAnnot_parsed", default=None, help="MatchAnnot parsed output (default: None)")
    parser.add_argument("-m2", "--sqanti_class", default=None, help="SQANTI classification file (default: None)")
    parser.add_argument("--demux", default=None, help="Demuxed FL count file - if provided, will be used instead of the <input_prefix>.abundance.txt file")
    parser.add_argument("--include_single_exons", default=False, action="store_true", help="Include single exons (default: OFF)")

    args = parser.parse_args()
    make_file_for_subsample(args.input_prefix, args.output_prefix, args.demux, args.matchAnnot_parsed, args.sqanti_class, args.include_single_exons)



