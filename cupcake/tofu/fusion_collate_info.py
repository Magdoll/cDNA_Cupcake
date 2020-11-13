#!/usr/bin/env python

import os, re, sys
from collections import defaultdict
from csv import DictReader, DictWriter
from Bio import SeqIO
from cupcake.io.GFF import collapseGFFReader

fusion_pbid = re.compile('PBfusion.(\d+).(\d+)')
"""
Run after fusion_finder.py + SQANTI3 classification
"""

FIELDS = ['UniqueID', 'FusionName', 'LeftGeneName', 'LeftGeneID', 'LeftBreakpoint', 'LeftFlankingSequence',
          'RightGeneName', 'RightGeneID', 'RightBreakpoint', 'RightFlankingSequence',
          'JunctionSupport', 'SpanningReads', 'ReadCountScore',
          'Sequence', 'LeftORF', 'RightORF', 'LeftExonCount', 'RightExonCount',
          'LeftCDSExonCount', 'RightCDSExonCount']


def get_breakpoint_n_seq(r1, r2, genome_dict=None, flanking_size=50):
    if r1.strand == '+':
        left_breakpoint = "{0}:{1}:+".format(r1.chr, r1.end)
        if genome_dict is not None:
            left_seq = str(genome_dict[r1.chr][r1.end-flanking_size:r1.end].seq)
        else:
            left_seq = 'NA'
    else:
        left_breakpoint = "{0}:{1}:-".format(r1.chr, r1.start)
        if genome_dict is not None:
            left_seq = str(genome_dict[r1.chr][r1.start:r1.start+flanking_size].reverse_complement().seq)
        else:
            left_seq = 'NA'
    if r2.strand == '+':
        right_breakpoint = "{0}:{1}:+".format(r2.chr, r2.start)
        if genome_dict is not None:
            right_seq = str(genome_dict[r2.chr][r2.start:r2.start+flanking_size].seq)
        else:
            right_seq = 'NA'
    else:
        right_breakpoint = "{0}:{1}:-".format(r2.chr, r2.end)
        if genome_dict is not None:
            right_seq = str(genome_dict[r2.chr][r2.end-flanking_size:r2.end].reverse_complement().seq)
        else:
            right_seq = 'NA'
    return left_breakpoint, left_seq, right_breakpoint, right_seq

def collate_info(fusion_prefix, class_filename, genepred_filename,
                 total_fl_count=None,
                 config_filename=None,
                 genome_dict=None,
                 cds_gff_filename=None,
                 min_fl_count=2):

    global_info = {}   # holding information for general information
    if config_filename is not None:
        print("Reading config file {0}...".format(config_filename), file=sys.stdout)
        for line in open(config_filename):
            k, v = line.strip().split('=')
            global_info[k] = v

    gene_to_id = {} # gene name --> ensembl ID
    for line in open(genepred_filename):
        raw = line.strip().split()
        gene_to_id[raw[11]] = raw[0]


    d = defaultdict(lambda: {}) # PBfusion.X --> isoform index -> sqanti3 record
    orf_dict = {}
    # read SQANTI3 classification file
    for r in DictReader(open(class_filename), delimiter='\t'):
        m = fusion_pbid.match(r['isoform'])
        if m is None:
            print("ERROR: fusion pbid must follow format `PBfusion.X.Y`. Abort!", file=sys.stderr)
            sys.exit(-1)
        gene_index, isoform_index = m.group(1), m.group(2)
        d[gene_index][isoform_index] = r
        orf_dict[r['isoform']] = r['ORF_seq']

    # get sequences
    seq_dict = dict((r.id.split('|')[0], r.seq) for r in SeqIO.parse(open(fusion_prefix + '.rep.fa'),'fasta'))

    # get count information
    count_d = defaultdict(lambda: 'NA')
    count_filename = fusion_prefix + '.abundance.txt'
    if os.path.exists(count_filename):
        for r in DictReader(open(count_filename), delimiter='\t'):
            count_d[r['pbid']] = int(r['count_fl'])

    if total_fl_count is None:
        print("Total FL count not given --- using the sum FL count from fusions only instead.", file=sys.stdout)
        total_fl_count = sum(count_d.values())

    # get breakpoint information
    gff_d = defaultdict(lambda: {}) # PBfusion.X --> isoform index -> sqanti3 record
    if cds_gff_filename is None:
        gff_filename = fusion_prefix + '.gff'
    else:
        gff_filename = cds_gff_filename

    for r in collapseGFFReader(gff_filename):
        m = fusion_pbid.match(r.seqid)
        if m is None:
            print("ERROR: fusion pbid in {0} must follow format `PBfusion.X.Y`. Abort!".format(gff_filename), file=sys.stderr)
            sys.exit(-1)
        gene_index, isoform_index = m.group(1), int(m.group(2))
        gff_d[gene_index][isoform_index] = r
        if r.strand not in ('+','-'):
            print("ERROR: fusion {0} did not specify strand in {1}! Abort!".format(r.seqid, gff_filename))
            sys.exit(-1)

    fields2 = list(global_info.keys()) + FIELDS
    f = open(fusion_prefix + '.annotated.txt', 'w')
    f_bad = open(fusion_prefix + '.annotated_ignored.txt', 'w')
    writer = DictWriter(f, fields2, delimiter=',')
    writer.writeheader()
    writer_bad = DictWriter(f_bad, fields2, delimiter=',')
    writer_bad.writeheader()

    for gene_index, iso_dict in d.items():
        iso_dict = list(iso_dict.items())  # (isoform index, classification record)
        iso_dict.sort(key=lambda x: x[0])
        has_novel = any(r['associated_gene'].startswith('novelGene') or r['associated_gene']=='' for junk,r in iso_dict)
        pbid = 'PBfusion.' + str(gene_index)

        gff_info = list(gff_d[gene_index].items())
        gff_info.sort(key=lambda x: x[0])

        rec1 = gff_info[0][1]
        rec2 = gff_info[-1][1]
        left_breakpoint, left_seq, right_breakpoint, right_seq = \
            get_breakpoint_n_seq(rec1, rec2, genome_dict)
        left_exon_count = len(rec1.ref_exons)
        right_exon_count = len(rec2.ref_exons)
        gene1 = iso_dict[0][1]['associated_gene']
        gene2 = iso_dict[-1][1]['associated_gene']

        if cds_gff_filename is not None:
            left_cds_exon_count = len(rec1.cds_exons)
            right_cds_exon_count = len(rec2.cds_exons)
        else:
            left_cds_exon_count = 'NA'
            right_cds_exon_count = 'NA'

        left_orf, right_orf = 'NA', 'NA'
        if orf_dict is not None:
            seqid1 = gff_info[0][1].seqid
            seqid2 = gff_info[-1][1].seqid
            left_orf = orf_dict[seqid1]
            right_orf = orf_dict[seqid2]

        info = {'UniqueID': pbid,
                'FusionName': "--".join([_r['associated_gene'] for (_index,_r) in iso_dict]),
                'LeftGeneName': gene1,
                'LeftGeneID': gene_to_id[gene1] if gene1 in gene_to_id else 'NA',
                'LeftBreakpoint': left_breakpoint,
                'LeftFlankingSequence': left_seq,
                'RightGeneName': gene2,
                'RightGeneID': gene_to_id[gene2] if gene2 in gene_to_id else 'NA',
                'RightBreakpoint': right_breakpoint,
                'RightFlankingSequence': right_seq,
                'JunctionSupport': 'NA',
                'SpanningReads': count_d[pbid],
                'ReadCountScore': count_d[pbid]*(10**6)/total_fl_count  if count_d[pbid] is not 'NA' else 'NA',
                'Sequence': seq_dict[pbid],
                'LeftORF': left_orf,
                'RightORF': right_orf,
                'LeftExonCount': left_exon_count,
                'RightExonCount': right_exon_count,
                'LeftCDSExonCount': left_cds_exon_count,
                'RightCDSExonCount': right_cds_exon_count}
        info.update(global_info)
        if has_novel or \
                gene1==gene2 or \
                (info['SpanningReads']!='NA' and info['SpanningReads'] < min_fl_count):
            writer_bad.writerow(info)
        else:
            writer.writerow(info)

    f.close()

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("fusion_prefix", help="Prefix for fusion files (ex: my.fusion)")
    parser.add_argument("class_filename", help="SQANTI3 classification filename")
    parser.add_argument("genepred_filename", help="GenePred filename used by SQANTI3 classification")
    parser.add_argument("--total_fl_count", type=int, default=None, help="(optional) Total FL count used to normalize fusion counts")
    parser.add_argument("--config", help="(optional) Additional information to include in the output")
    parser.add_argument("--genome", help="(optional) Reference genome")
    parser.add_argument("--min_fl_count", type=int, default=2, help="Minimum FL count (default: 2)")

    args = parser.parse_args()

    if args.genome is not None:
        genome_dict = SeqIO.to_dict(SeqIO.parse(open(args.genome), 'fasta'))
        print("Finished reading reference genome {0}.".format(args.genome))
    else:
        genome_dict = None

    collate_info(args.fusion_prefix, args.class_filename, args.genepred_filename,
                 total_fl_count=args.total_fl_count,
                 config_filename=args.config,
                 genome_dict=genome_dict,
                 min_fl_count=args.min_fl_count)