#!/usr/bin/env python
import os, sys
from cupcake.io.GFF import write_collapseGFF_format
from cupcake.io.BioReaders import GMAPSAMReader

def main():
    from argparse import ArgumentParser

    parser = ArgumentParser("Convert SAM to collapsed GFF format")
    parser.add_argument("sam_filename")

    args = parser.parse_args()

    if not args.sam_filename.endswith('.sam'):
        print >> sys.stderr, "Only accepts files ending in .sam. Abort!"
        sys.exit(-1)

    prefix = args.sam_filename[:-4]
    output_gff = prefix + '.collapsed.gff'

    with open(output_gff, 'w') as f:
        reader = GMAPSAMReader(args.sam_filename, True)
        for r in reader:
            if r.sID == '*': continue 
            r.strand = r.flag.strand
            r.seqid = r.qID
            r.chr = r.sID
            r.ref_exons = r.segments
            r.start = r.sStart
            r.end = r.sEnd
            r.cds_exons = None
            write_collapseGFF_format(f, r)

    print >> sys.stderr, "Output written to {0}.".format(output_gff)

if __name__ == "__main__":
    main()
