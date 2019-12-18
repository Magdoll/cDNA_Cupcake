#!/usr/bin/env python
__author__ = 'etseng@pacb.com'

import os, sys, re
from bx.intervals.cluster import ClusterTree
import cupcake.io.GFF as GFF
from Bio import SeqIO

"""
Given a GFF file and selected loci (ex: PB.45),
create a fake genome that is the concatenation of all seen exons.
"""

extra_bp_around_junctions = 50 # get this much around junctions to be safe AND to not screw up GMAP who doesn't like microintrons....
__padding_before_after__ = 10 # get this much before and after the start

def make_fake_genome(genome_filename, gff_filename, ref_chr, ref_start, ref_end, ref_strand, output_prefix, output_name, genome_d=None):
    if genome_d is None:
        print("Reading genome file {0}...".format(genome_filename), file=sys.stderr)
        d = SeqIO.to_dict(SeqIO.parse(open(genome_filename),'fasta'))
    else:
        d = genome_d

    print("Reading GFF file {0}...".format(gff_filename), file=sys.stderr)
    good = []
    reader = GFF.collapseGFFReader(gff_filename)
    for r in reader:
        if r.chr==ref_chr and r.strand==ref_strand and \
                (ref_start <= r.start < r.end <= ref_end) \
            and len(r.ref_exons) > 1:
            print("Adding {0} to fake genome.".format(r.seqid), file=sys.stderr)
            good.append(r)

    if len(good) == 0:
        print("Did not find any transcripts strictly within {0}:{1}-{2} on strand {3}. Abort!".format(\
            ref_chr, ref_start, ref_end, ref_strand), file=sys.stderr)
        sys.exit(-1)

    c = ClusterTree(0, 0)
    for r in good:
        for e in r.ref_exons:
            c.insert(e.start-extra_bp_around_junctions, e.end+extra_bp_around_junctions, 1)

    regions = [(a,b) for (a,b,junk) in c.getregions()]
    regions[0] = (regions[0][0]-__padding_before_after__, regions[0][1])
    regions[-1] = (regions[-1][0], regions[-1][1]+__padding_before_after__)

    with open(output_prefix+'.fasta', 'w') as f:
        f.write(">" + output_name + "\n")
        for a,b in regions:
            f.write(str(d[r.chr][a:b].seq))
        f.write("\n")
        f.close()

    # for mapping, write <0-based index on fake genome>, <ref chrom>, <0-based index on ref genome>
    with open(output_prefix+'.mapping.txt', 'w') as f:
        i = 0
        for a,b in regions:
            for j in range(a, b):
                f.write("{0},{1},{2}\n".format(i, ref_chr, j))
                i += 1

        with open(output_prefix+'.pbids.txt', 'w') as f:
            f.write("\n".join(r.seqid for r in good)+'\n')

    print("Output written to {0}.fasta, {0}.mapping.txt, {0}.pbids.txt.".format(output_prefix), file=sys.stderr)


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("genome_filename")
    parser.add_argument("gff_filename")
    parser.add_argument("--locus", required=True, help="locus in format <chr>:<start>-<end>")
    parser.add_argument("--strand", required=True, choices=['+', '-'], help="strand of locus")
    parser.add_argument("-o", "--output_prefix", help="Output prefix")

    args = parser.parse_args()

    rex = re.compile('(\S+):(\d+)-(\d+)')
    m = rex.match(args.locus)
    if m is None:
        print("{0} is not a defined chr location! Abort.".format(args.locus), file=sys.stderr)
        sys.exit(-1)

    ref_chr = m.group(1)
    ref_start = int(m.group(2)) - 1 # make it 0-based
    ref_end = int(m.group(3)) # keep it 1-based

    if not os.path.exists(args.genome_filename):
        print("Genome file {0} does not exist! Abort.".format(args.genome_filename), file=sys.stderr)
        sys.exit(-1)
    if not os.path.exists(args.gff_filename):
        print("GFF {0} does not exist! Abort.".format(args.gff_filename), file=sys.stderr)
        sys.exit(-1)

    make_fake_genome(args.genome_filename, args.gff_filename, ref_chr, ref_start, ref_end, args.strand, args.output_prefix)
