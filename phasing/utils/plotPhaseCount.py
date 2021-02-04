# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 13:14:28 2020

@author: derek.bickhart-adm
"""

import matplotlib
from matplotlib import pyplot as plt
matplotlib.use('Agg')
from matplotlib.collections import BrokenBarHCollection
from matplotlib import cm
from itertools import cycle
from collections import defaultdict
import argparse
import pandas
import numpy as np
import pysam

def arg_parse():
    parser = argparse.ArgumentParser(
            description = "A tool to plot bin and contig level read depth differences in strain assignment"
            )
    parser.add_argument('-f', '--fai', 
                        help="Input reference fasta index file for the bin",
                        required=True, type=str
                        )
    parser.add_argument('-o', '--output',
                        help="Output file basename. Output files are {output}.wins and {output}.pdf",
                        required=True, type=str,
                        )
    parser.add_argument('-b', '--bam', 
                        help="Input CCS read depth bam file",
                        required=True, type=str
                        )
    parser.add_argument('-h', '--human', 
                        help="Input human-readable variant call file",
                        required=True, type=str
                        )
    parser.add_argument('-i', '--binsize',
                        help="Bin size in bases [5000 bp]",
                        type = int, default=5000
                        )
    return parser.parse_args(), parser
    
def main(args, parser):
    # Get the contig length list
    ctglens = dict()
    with open(args.fai, 'r') as fai:
        for l in fai:
            s = l.rstrip().split()
            ctglens[s[0]] = s[1]
            
    # Create windows 
    winlist = defaultdict(list)
    # offset bp to add for stitching contigs together in one line
    ctgoffset = dict()
    lastbp = 0
    for c in ctglens:
        ctgoffset[c] = lastbp + 100
        for i in range(0, ctglens[c], args.binsize):
            winlist[c].append(window(c, i, i + args.binsize))
        lastbp += ctglens[c]
        
    # read each sam region and count the reads
    with pysam.AlignmentFile(args.bam, 'rb') as bamfile:
        for c, w in winlist.items():
            for i, win in enumerate(w):
                count = 0
                for s in bamfile.fetch(c, win.start, win.end):
                    if s.is_secondary:
                        continue
                    count += 1
                winlist = updateWin(winlist, c, i, count)
                
    # Now, read in the human readable text file and process that 
    hapset = set()
    with open(args.human, 'r') as human:
        human.readline()
        for l in human:
            s = l.rstrip().split()
            # determine where the contig start falls
            for i, win in enumerate(winlist[s[2]]):
                if int(s[3]) < win.end and int(s[3]) >= win.start:
                    winlist = updateWin(winlist, s[2], i, int(s[6]), s[4])
                    print(f'Updating window: {s[2]} {win.start} {win.end} to {s[6]} for Hap {s[4]}')
                    hapset.add(s[4])
                    
    # OK, data is in! Let's try plotting
    raw = defaultdict(list)
    bars = list()
    for c, w in winlist.items():
        bars.append([ctgoffset[c], ctglens[c]])
        for win in winlist:
            for h in hapset:
                raw["contig"].append(c)
                raw["start"].append(win.start + ctgoffset[c])
                raw["end"].append(win.end + ctgoffset[c])
                raw["hap"].append(h)
                raw["count"].append(win.getCount(h))
                
    df = pandas.DataFrame(raw)
    df.to_csv(args.output + '.wins', sep='\t', header=True)
    
    fig = plt.figure(figsize=(6,8))
    ax = df[['start', 'hap', 'count']].plot.area(x='start', y='count', colormap='viridis')
    
    ax.add_collection(BrokenBarHCollection(bars, [-1, 1], facecolors=plt.get_cmap('tab20')))
    ax.axis('tight')
    plt.savefig(args.output + '.pdf')
    
    

            
            
    
def updateWin(winlist, contig, winidx, count, haplotype = 'REF'):
    winlist[contig].count[haplotype] = count
    return winlist
    

class window:
    
    def __init__(self, contig, start, end):
        self.contig = contig
        self.start = start 
        self.end = end 
        self.count = defaultdict(int)
        
    def getCount(self, hap):
        if hap in self.count:
            return self.count[hap]
        else:
            return 0
    

if __name__ == "__main__":
    args, parser = arg_parse()
    main(args, parser)
