#!/usr/bin/env python3
# Roger Volden

'''
Meant to plot the number of reads per cell barcode,
can also plot umis per cell

Usage:
    python3 plot_knees.py \
            -b skera.fltnc.bam
            -o plots/skera
            -t {cbc,umi,both}
'''

import argparse
import pysam
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--bam', '-b', type=str, required=True, help='Input bam'
    )
    parser.add_argument(
        '--output', '-o', type=str, required=True, help='Output png prefix'
    )
    parser.add_argument(
        '--type', '-t', type=str,
        choices=['cbc', 'umi', 'both'], default='both',
        help='Type of plot to output (default both)'
    )
    parser.add_argument(
        '--max_cells', '-m', default=-1, type=int,
        help='Force an x axis maximum instead of the mpl default'
    )
    return parser.parse_args()

def read_counts(counts):
    '''
    Reads in the dedup csv file and keeps track
    of how many reads per cell barcode.

    read_counts = {'cell_bc_1': [read_count, {'molecule/1', ...}], ...}
    '''
    read_counts = {}
    pysam.set_verbosity(0)
    bam = pysam.AlignmentFile(counts, 'rb', check_sq=False)
    for line in tqdm(bam.fetch(until_eof=True), desc='Reading bam'):
        line = line.tostring(bam)
        line = line.rstrip().split('\t')
        read_id = line[0]
        cbc, umi = '', ''
        for tag in line[11:]:
            if tag.startswith('XC:Z:'):
                cbc = tag[5:]
            if tag.startswith('XM:Z:'):
                umi = tag[5:]
        if not cbc or not umi:
            continue
        if cbc not in read_counts:
            read_counts[cbc] = [1, set([umi])]
        else:
            read_counts[cbc][1].add(umi)
            read_counts[cbc][0] += 1
    bam.close()
    return read_counts

def plot_rpc(counts, max_cells, output):
    plt.figure(figsize=(5, 3))
    plt.style.use('clean.mplstyle')
    c = plt.axes([0.125, 0.125, 0.8, 0.8])

    # only look at the read counts per cell
    y = sorted([x[0] for x in list(counts.values())], reverse=True)
    print(f'Total reads up to cutoff: {sum(y[:max_cells])}')
    if max_cells < 0:
        c.plot(range(len(y)), np.log2(y), lw=1, color='black')
    else:
        c.plot(
            range(len(y[:max_cells])), 
            np.log2(y[:max_cells]), 
            lw=1, color='black'
        )

    c.set_xlabel('Cell #')
    c.set_ylabel(r'log$_2$(# of reads)')
    c.set_title('Reads per cell')

    output += '.rpc.png'
    plt.savefig(output, dpi=600)

def plot_upc(counts, max_cells, output):
    plt.figure(figsize=(5, 3))
    plt.style.use('clean.mplstyle')
    c = plt.axes([0.125, 0.125, 0.8, 0.8])

    # look at the umi set
    y = sorted([len(x[1]) for x in list(counts.values())], reverse=True)
    print(f'Total UMIs up to cutoff: {sum(y[:max_cells])}')
    if max_cells < 0:
        c.plot(range(len(y)), np.log2(y), lw=1, color='black')
    else:
        c.plot(
            range(len(y[:max_cells])), 
            np.log2(y[:max_cells]), 
            lw=1, color='black'
        )

    c.set_xlabel('Cell #')
    c.set_ylabel(r'log$_2$(# of reads)')
    c.set_title('UMIs per cell')

    output += '.upc.png'
    plt.savefig(output, dpi=600)

def main(args):
    counts = read_counts(args.bam)
    if args.type == 'cbc':
        plot_rpc(counts, args.max_cells, args.output)
    elif args.type == 'umi':
        plot_upc(counts, args.max_cells, args.output)
    else:
        plot_rpc(counts, args.max_cells, args.output)
        plot_upc(counts, args.max_cells, args.output)

if __name__ == '__main__':
    args = parse_args()
    main(args)
