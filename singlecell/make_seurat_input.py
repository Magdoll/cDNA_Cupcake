#!/usr/bin/env python3
# Roger Volden

'''
Usage:
    # to make isoform inputs
    python3 make_seurat_input.py \
        -i dedup.annotated.csv.gz \
        -t isoform \
        -o path/to/output

    # to make gene inputs
    python3 make_seurat_input.py \
        -i dedup.annotated.csv.gz \
        -t gene \
        -a gencode.v29.gtf \
        -o path/to/output

Output structure:
    path/to/output
    ├── genes_seurat
    │   ├── barcodes.tsv
    │   ├── genes.tsv
    │   └── matrix.mtx
    └── isoforms_seurat
        ├── barcodes.tsv
        ├── genes.tsv
        └── matrix.mtx
'''

import argparse
import gzip
import os
import sys
import shutil

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--input', '-i', type=str,
        help='Input dedup annotated csv file (can be gzipped)'
    )
    parser.add_argument(
        '--type', '-t', type=str,
        choices=['gene', 'isoform'], default='isoform',
        help='Produce output for gene expression or isoform expression (default: isoform)'
    )
    parser.add_argument(
        '--annotation', '-a', type=str,
        help='GTF file that has the gene or isoform information. If doing both gene and \
              isoform information, put the gencode GTF first (eg. gencode.gtf,dedup.anno.gtf).'
    )
    parser.add_argument(
        '--output', '-o', type=str, default=os.getcwd(),
        help='Output path. Will make genes/isoforms_seurat in this dir'
    )
    parser.add_argument(
        '--keep_RiboMito', '-r', action='store_true', default=False,
        help='Use to keep ribosomal and mitochondrial genes in the output (default: False)'
    )
    parser.add_argument(
        '--keep_novel', '-n', action='store_true', default=False,
        help='Use to keep novel genes in the output (default: False)'
    )
    parser.add_argument(
        '--clobber', '-c',
        action='store_true', default=False,
        help='Use to destroy existing output (default: False)'
    )
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    return parser.parse_args()

def check_args(args):
    if args.type == 'gene' and not args.annotation:
        print('You need to provide a GTF/GFF file', file=sys.stderr)
        exit(1)
    if not args.input:
        print('You need to provide an input csv file', file=sys.stderr)
        exit(1)

def read_gtf(gtf, args):
    '''
    Reads in GTF/GFF files.
    Handles the normal case where you have a gencode annotation.
    This will make a dictionary of gene ID to gene name and write
    it out to the genes.tsv file.
    '''
    out_path = args.output
    id_dict, gene_dict = {}, {}
    with open(gtf) as f:
        for line in f:
            if line[0] == '#':
                continue
            line = line.rstrip().split('\t')
            line = [x.strip() for x in line[8].split(';')]
            id, name = '', ''
            for section in line:
                if section.startswith('gene_id'):
                    id = section.split('"')[1]
                if section.startswith('gene_name'):
                    name = section.split('"')[1]
                if id and name:
                    break
            id_dict[id] = name

    if not os.path.isdir(out_path):
        print(f"Couldn't find {out_path}, making it now", file=sys.stderr)
        os.mkdir(out_path)
    out_path += 'genes_seurat/'
    if os.path.isdir(out_path) and os.listdir(out_path) and not args.clobber:
        print(f'A previous output exists in {out_path}. Use --clobber (-c) to overwrite', file=sys.stderr)
        exit(1)
    elif os.path.isdir(out_path) and os.listdir(out_path) and args.clobber:
        shutil.rmtree(out_path)

    os.mkdir(out_path)

    genes_out = open(out_path + 'genes.tsv', 'w+')
    index = 1
    print(f'Writing genes to {out_path}genes.tsv', file=sys.stderr)
    for gene_id, gene_name in id_dict.items():
        print(gene_id + '\t' + gene_name, file=genes_out)
        gene_dict[gene_name] = index
        index += 1
    genes_out.close()
    return gene_dict

def read_csv(csv, gene_dict, args):
    '''
    Reads in the csv file to get the list of cell barcodes and
    the per cell gene or isoform information.
    The per cell info is going to depend on the exp_type.
    '''
    out_path = args.output
    cell_bcs, zipped, first = {}, False, True
    if args.input.endswith('.gz'):
        csv_fh = gzip.open(args.input, 'rb')
        zipped = True
    else:
        csv_fh = open(args.input, 'r')

    pb = True if args.type == 'isoform' else False

    count_dict = {} # cell_bc: {gene_id: 1, pb_id: 1}
    index = 1
    gene_idx = 1
    total_counts = 0
    for line in csv_fh:
        if zipped:
            line = line.decode()
        line = line.rstrip().split('\t')
        if first:
            first = False
            continue
        # keep track of the cell barcodes and their order
        cell_bc = line[10]
        if cell_bc not in cell_bcs:
            cell_bcs[cell_bc] = index
            index += 1
        # gene_id might be 'novel', also could be ensembl
        pb_iso, gene_id, gene_name = line[1], line[3], line[4]

        if not args.keep_RiboMito:
            if gene_name.startswith('MT-') \
                    or gene_name.startswith('RPL') \
                    or gene_name.startswith('RPS'):
                continue
        if not args.keep_novel and gene_id == 'novel':
            continue

        if pb:
            if pb_iso not in gene_dict:
                gene_dict[pb_iso] = (gene_name, gene_idx)
                gene_idx += 1
        if cell_bc not in count_dict:
            count_dict[cell_bc] = {pb_iso: 0}#, gene_id: 0}
        if pb and pb_iso not in count_dict[cell_bc]:
            count_dict[cell_bc][pb_iso] = 0
        if not pb and gene_name not in count_dict[cell_bc]:
            count_dict[cell_bc][gene_name] = 0
        if pb:
            count_dict[cell_bc][pb_iso] += 1
        else:
            count_dict[cell_bc][gene_name] += 1

        total_counts += 1

    if not os.path.isdir(out_path):
        print(f"Couldn't find {out_path}, making it now", file=sys.stderr)
        os.mkdir(out_path)
    out_path += args.type + 's_seurat/'
    if pb:
        if os.path.isdir(out_path) and os.listdir(out_path) and not args.clobber:
            print(f'A previous output exists in {out_path}. Use --clobber (-c) to overwrite', file=sys.stderr)
            exit(1)
        elif os.path.isdir(out_path) and os.listdir(out_path) and args.clobber:
            shutil.rmtree(out_path)

    if not os.path.isdir(out_path):
        os.mkdir(out_path)

    if pb:
        gene_file = args.output + args.type + 's_seurat/genes.tsv'
        print(f'Writing genes to {gene_file}', file=sys.stderr)
        gene_fh = open(gene_file, 'w+')
        for pbid, tup in gene_dict.items():
            print(f'{pbid}\t{tup[0]}', file=gene_fh)
        gene_fh.close()

    bc_file = args.output + args.type + 's_seurat/barcodes.tsv'
    print(f'Writing barcodes to {bc_file}', file=sys.stderr)
    bc_fh = open(bc_file, 'w+')
    for bc in sorted(cell_bcs.items(), key=lambda x:x[1]):
        print(f'{bc[0]}-1', file=bc_fh)
    bc_fh.close()

    num_cells, num_genes = len(cell_bcs), len(gene_dict)
    mtx_file = args.output + args.type + 's_seurat/matrix.mtx'
    print(f'Writing matrix to {mtx_file}', file=sys.stderr)
    mtx_fh = open(mtx_file, 'w+')
    print('%%MatrixMarket matrix coordinate real general\n%', file=mtx_fh)
    print(f'{num_genes} {num_cells} {total_counts}', file=mtx_fh)
    for cell, counts in count_dict.items():
        for gene, count in counts.items():
            if gene not in gene_dict:
                continue
            if pb:
                gene_idx, cell_idx = gene_dict[gene][1], cell_bcs[cell]
            else:
                gene_idx, cell_idx = gene_dict[gene], cell_bcs[cell]
            print(f'{gene_idx} {cell_idx} {count}', file=mtx_fh)
    mtx_fh.close()
    csv_fh.close()

def main(args):
    if not args.output.endswith('/'):
        args.output += '/'

    pb = True if args.type == 'isoform' else False
    if not pb:
        gene_dict = read_gtf(args.annotation, args)
    else:
        gene_dict = {}
    read_csv(args.input, gene_dict, args)

if __name__ == '__main__':
    args = parse_args()
    main(args)
