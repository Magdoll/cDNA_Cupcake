#!/usr/bin/env python3

import os, sys
from csv import DictReader
import pysam

def paint_bam_post_phaser(input_bam, output_bam, read_info, chrom, start, end):
    reader = pysam.AlignmentFile(input_bam, 'rb')
    fout = pysam.AlignmentFile(output_bam, 'wb', header=reader.header)
    for r in reader.fetch(chrom, start, end):
        d = r.to_dict()
        newtags = []
        for k in d['tags']:
            if not k.startswith('RG:Z'):
                newtags.append(k)
        if r.qname not in read_info:
            newtags.append('RG:Z:unassigned')
        else:
            hapstr = read_info[r.qname]
            newtags.append('RG:Z:' + hapstr)
        d['tags'] = newtags
        fout.write(pysam.AlignedSegment.from_dict(d, r.header))
    fout.close()

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("input_bam")
    parser.add_argument("output_bam")
    parser.add_argument("read_hap_info", help="Human readable read-to-hap info file from Isophase/Magphase output")
    parser.add_argument("-c", "--chrom", required=True, help='Chromosome')
    parser.add_argument("-s", "--start", required=True, type=int, help="Start location")
    parser.add_argument("-e", "--end", required=True, type=int, help="End location")

    args = parser.parse_args()

    read_info = {}
    for r in DictReader(open(args.read_hap_info), delimiter='\t'):
        read_info[r['read_id']] = r['haplotype']

    paint_bam_post_phaser(args.input_bam, args.output_bam, read_info, args.chrom, args.start, args.end)