#!/usr/bin/env python
import os, sys
from collections import defaultdict, Counter
from cupcake.io.GFF import write_collapseGFF_format
from cupcake.io.BioReaders import GMAPSAMReader

def main():
    from argparse import ArgumentParser

    parser = ArgumentParser("Convert SAM to collapsed GFF format")
    parser.add_argument("sam_filename")

    args = parser.parse_args()

    if not args.sam_filename.endswith('.sam'):
        print("Only accepts files ending in .sam. Abort!", file=sys.stderr)
        sys.exit(-1)

    prefix = args.sam_filename[:-4]
    output_gff = prefix + '.collapsed.gff'

    ids_seen = Counter()
    with open(output_gff, 'w') as f:
        reader = GMAPSAMReader(args.sam_filename, True)
        for r in reader:
            if r.sID == '*': continue 
            invalid_flag = False
            # check if any of the exons are invalid
            for (s,e) in r.ref_exons:
                if s >= e:
                    print("{0} has an invalid exon {1}:{2}-{3}. This record will not be output in GFF!".format(r.qID, r.sID, s, e))
                    invalid_flag = True
                    break
            if invalid_flag: continue
            
            if r.qID in ids_seen:
                i = ids_seen[r.qID] + 1
                print("{0} has multiple alignments! Adding a suffix {1} to {2}:{3}-{4}".format(r.qID, i, r.sID, r.sStart, r.sEnd))
                ids_seen[r.qID] += 1
                r.qID += ".dup" + str(i)
            else:
                ids_seen[r.qID] += 1
            r.strand = r.flag.strand
            r.geneid = r.qID
            r.seqid = r.qID
            r.chr = r.sID
#            r.ref_exons = r.segments
            r.start = r.sStart
            r.end = r.sEnd
            r.cds_exons = None
            write_collapseGFF_format(f, r)

    print("Output written to {0}.".format(output_gff), file=sys.stderr)

if __name__ == "__main__":
    main()
