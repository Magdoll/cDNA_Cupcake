#!/usr/bin/env python
import os, sys
from csv import DictReader, DictWriter
from collections import Counter, defaultdict
from Bio.Seq import Seq

def read_dropseq_clean_report(report_filename):
    """
    EXAMPLE
    # FILTER_AMBIGUOUS=true
    # MIN_UMIS_PER_CELL=20
    # UMI_BIAS_THRESHOLD=20
    # EDIT_DISTANCE=20
    #
    # TOTAL_BARCODES_TESTED=996206
    # BARCODES_COLLAPSED=12546
    # ESTIMATED_UMIS_COLLAPSED=1107912
    # AMBIGUOUS_BARCODES=8206
    # ESTIMATED_AMBIGUOUS_UMIS=427356
    # POLY_T_BIASED_BARCODES=3374
    # POLY_T_BIASED_BARRCODE_UMIS=568661
    # POLY_T_POSITION=8
    intended_barcode        neighbor_barcode        intended_size   neighbor_size   position        intended_base   neighbor_base   repaired
    GGTTGGTCGCGG    GGTTGGTCTCGG    5458    20      9       G       T       FALSE
    TGGAGCTGGCGC    TGGCGCTGGCGC    3826    190     4       A       C       TRUE
    GGCGGCACGGCC    GACGGCACGGCC    16997   69      2       G       A       FALSE
    GAAGCCAGAGGT    GAAGCCAGAAGT    7037    47      10      G       A       FALSE
    GTGACGGACGGT    GTGACCGACGGT    8959    69      6       G       C       FALSE
    """
    f = open(report_filename)
    while True:
        cur_pos = f.tell()
        if not f.readline().startswith('#'): break
    f.seek(cur_pos)

    bc_repair_dict = {}
    reader = DictReader(f, delimiter='\t')
    for r in reader:
        if r['repaired']=='TRUE':
            seq_from = Seq(r['neighbor_barcode']).reverse_complement()
            seq_to   = Seq(r['intended_barcode']).reverse_complement()
            bc_repair_dict[seq_from] = seq_to
    return bc_repair_dict

def read_dropseq_synthesis_report(report_filename, bc_repair_dict=None):
    """
    EXAMPLE
    intended_sequence       related_sequences       num_related     deleted_base    deleted_base_pos        non_incorporated_rate   intended
_UMIs   related_median_UMIs     intended_TBias  related_median_TBias
    NA      CGCAGCTCTGAG    1       NA      NA      NA      NA      20      NA      1
    NA      GCTGGACTTACA    1       NA      NA      NA      NA      22      NA      0.95
    CTCCACTGGAAA    CTCCACTGAAAA:CTCCACTGAAAC:CTCCACTGAAAG:CTCCACTGAAAT     4       G       9       0.47    4851    1039    0.3     0.98
    TTTAATATGGAT    TTAATATGGATG:TTAATATGGATT
    """
    f = open(report_filename)
    while True:
        cur_pos = f.tell()
        if not f.readline().startswith('#'): break
    f.seek(cur_pos)

    if bc_repair_dict is None:
        bc_repair_dict = {}
    reader = DictReader(f, delimiter='\t')
    for r in reader:
        if r['intended_sequence']!='NA':
            seq_to   = Seq(r['intended_sequence']).reverse_complement()
            for s in r['related_sequences'].split(':'):
                seq_from = Seq(s).reverse_complement()
                bc_repair_dict[seq_from] = seq_to
    return bc_repair_dict

def edit_distance(seq1, seq2):
    assert len(seq1)==len(seq2)
    diff = 0
    for i in range(len(seq1)):
        diff += seq1[i]!=seq2[i]
    return diff

def error_correct_BC_or_UMI(records, key, threshold=1):
    """
    :param records: should be list of records all from the same gene!
    """
    assert key in ('BC', 'UMI')
    merge_map = {}
    bc_count = Counter()
    for r in records: bc_count[r[key]] += 1

    # most common BC, in decreasing order
    bc_sorted = [bc for bc,count in bc_count.most_common()]

    i = 0
    while i < len(bc_sorted)-1:
        j = i + 1
        while j < len(bc_sorted):
            if edit_distance(bc_sorted[i], bc_sorted[j]) <= threshold:
                #merge j into i
                merge_map[bc_sorted[j]] = bc_sorted[i]
                bc_sorted.pop(j)
            else:
                j += 1
        i += 1

    #if len(merge_map) > 0:
    #    print merge_map
     ##   raw_input()
    return merge_map


def main(csv_filename, output_filename, shortread_bc={}, only_top_ranked=False, bc_repair_dict=None):

    reader = DictReader(open(csv_filename), delimiter='\t')

    FIELDS = reader.fieldnames + ['BC_ed', 'UMI_ed', 'BC_match', 'BC_top_rank']
    f = open(output_filename, 'w')
    writer = DictWriter(f, FIELDS, delimiter='\t')
    writer.writeheader()

    recs_by_gene = defaultdict(lambda: [])
    for r in reader:
        recs_by_gene[r['gene']].append(r)

    # error correct BCs by gene group
    for gene in recs_by_gene:
        recs_by_bc = defaultdict(lambda: [])

        if bc_repair_dict is not None: # has DropSeq BC cleaning report! Just use it!
            for r in recs_by_gene[gene]:
                if r['BC'] in bc_repair_dict:
                    r['BC_ed'] = bc_repair_dict[r['BC']]
                else:
                    r['BC_ed'] = r['BC']
                recs_by_bc[r['BC']].append(r)
        else:
            bc_merge_map = error_correct_BC_or_UMI(recs_by_gene[gene], 'BC')
            for r in recs_by_gene[gene]:
                if r['BC'] in bc_merge_map:
                    r['BC_ed'] = bc_merge_map[r['BC']]
                else:
                    r['BC_ed'] = r['BC']
                recs_by_bc[r['BC']].append(r)
        # now error correct by UMI
        for bc in recs_by_bc:
            umi_merge_map = error_correct_BC_or_UMI(recs_by_bc[bc], 'UMI')
            for r in recs_by_bc[bc]:
                if r['UMI'] in umi_merge_map:
                    r['UMI_ed'] = umi_merge_map[r['UMI']]
                else:
                    r['UMI_ed'] = r['UMI']

                BC_ed_rev = str(Seq(r['BC_ed']).reverse_complement())

                r['BC_match'] = 'Y' if BC_ed_rev in shortread_bc else 'N'
                r['BC_top_rank'] = 'Y' if (r['BC_match']=='Y' and shortread_bc[BC_ed_rev]=='Y') else 'N'

                if (not only_top_ranked or r['BC_top_rank']=='Y'):
                    writer.writerow(r)


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("input_csv", help="Input CSV")
    parser.add_argument("output_csv", help="Output CSV")
    parser.add_argument("--bc_rank_file", help="(Optional) cell barcode rank file from short read data")
    parser.add_argument("--only_top_ranked", action="store_true", default=False, help="(Optional) only output those that are top-ranked. Must have --bc_rank_file.")
    parser.add_argument("--dropseq_clean_report", help="Output from running DetectBeadSubstitutionErrors in DropSeq cookbook (ex: star_gene_exon_tagged_clean_substitution.bam_report.txt)")
    parser.add_argument("--dropseq_synthesis_report", help="Output from running DetectBeadSynthesisErrors in DropSeq cookbook (ex: star_gene_exon_tagged_clean_substitution_clean2.bam_report.txt)")

    args = parser.parse_args()

    shortread_bc = {}  # dict of cell barcode -> "Y" for top ranked
    if args.bc_rank_file is not None:
        reader = DictReader(open(args.bc_rank_file), delimiter='\t')
        for r in reader:
            shortread_bc[r['cell_barcode']] = r['top_ranked']
    else:
        if args.only_top_ranked:
            print("--bc_rank_file must be given if using --only_top_ranked!", file=sys.stderr)
            sys.exit(-1)

    bc_repair_dict = None
    if args.dropseq_clean_report is not None:
        bc_repair_dict = read_dropseq_clean_report(args.dropseq_clean_report)
    if args.dropseq_synthesis_report is not None:
        bc_repair_dict = read_dropseq_synthesis_report(args.dropseq_synthesis_report, bc_repair_dict)

    main(args.input_csv, args.output_csv, shortread_bc, args.only_top_ranked, bc_repair_dict)
