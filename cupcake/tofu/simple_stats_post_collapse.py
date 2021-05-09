#!/usr/bin/env python

import os, sys
from cupcake.io.GFF import collapseGFFReader

def main(input_prefix):
    input_gff = input_prefix + '.gff'
    if not os.path.exists(input_gff):
        print("Looking for input GFF {0} but not found! Abort!".format(input_gff))
        sys.exit(-1)

    f1 = open(input_prefix + '.simple_stats.txt', 'w')
    f1.write("pbid\tlocus\tlength\tgenomic_length\tnum_exon\n")
    f2 = open(input_prefix + '.exon_stats.txt', 'w')
    f2.write("pbid\texon_index\texon_size\tintron_size\n")
    for r in collapseGFFReader(input_gff):
        if len(r.ref_exons)==0: continue
        f1.write(r.seqid+'\t')
        if r.seqid.startswith('PB.'): f1.write(r.seqid.split('.')[1]+'\t')
        else: f1.write(r.geneid+'\t')
        sum_len = 0
        for i,e in enumerate(r.ref_exons):
            exon_len = e.end - e.start
            sum_len += exon_len
            f2.write("{0}\t{1}\t{2}\t".format(r.seqid, i+1, exon_len))
            if i == 0: f2.write("NA\n")
            else: 
                if len(r.ref_exons)>0: f2.write(str(e.start-r.ref_exons[i-1].end) + '\n')

        f1.write(str(sum_len)+'\t')
        f1.write(str(r.end-r.start)+'\t')
        f1.write(str(len(r.ref_exons))+'\n')

    f1.close()
    f2.close()
    print("Output written to: {0},{1}\n".format(f1.name, f2.name))

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("input_prefix", help="Input prefix, ex: hq.5merge.collapsed")

    args = parser.parse_args()
    main(args.input_prefix)
