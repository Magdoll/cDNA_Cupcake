__author__ = "etseng@pacb.com"
import pdb
"""
For Iso-Phase, selecting loci that has sufficient FL coverage to phase.

INPUT: 
   -- FLNC (fastq)
   -- Unique isoforms and genes (GFF), using PB.X.Y to denote gene vs isoforms relationship
   -- FLNC association to isoforms, (.read_stat.txt)
   -- Reference genome (fasta)
"""

import os, sys, re
from collections import defaultdict
from csv import DictReader
from Bio import SeqIO
from bx.intervals.cluster import ClusterTree
from cupcake.io.SeqReaders import LazyFastqReader, LazyFastaReader
from cupcake.io.GFF import collapseGFFReader

rex_flnc = re.compile('(m\S+_\d+_\d+\/\d+)\/ccs') # Sequel/Iso-Seq 3 format
rex_pbid = re.compile('(PB[fusion]?.\d+)[.\d+]?')  # this will handle both PB.X.Y and PBfusion.X


extra_bp_around_junctions = 50 # get this much around junctions to be safe AND to not screw up GMAP who doesn't like microintrons....
__padding_before_after__ = 10 # get this much before and after the start


def read_flnc_fafq(flnc_filename, fq_or_fa='fastq'):
    """
    Read FLNC fastq into a dict of zmw --> lazy file pointer
    """
    if fq_or_fa=='fastq':
        flnc_fafq_d = LazyFastqReader(flnc_filename)
    else:
        flnc_fafq_d = LazyFastaReader(flnc_filename)
    rich_zmws = set()

    for k in list(flnc_fafq_d.keys()):
        m = rex_flnc.match(k)
        if m is not None:
            zmw = m.group(1)
            flnc_fafq_d.d[zmw] = flnc_fafq_d.d[k]
        else:
            zmw = k
        rich_zmws.add(zmw)

    return flnc_fafq_d, rich_zmws


def read_read_stat(stat_filename, rich_zmws):
    """
    Read .read_stat.txt file
    :return: tally_by_loci -- dict  of {locus -- list of (isoform, zmw)}
             poor_zmws_not_in_rich --- list of ZMWs that were in .read_stat but missing in FLNC FASTQ. ONly happens if used diff CCS to get FASTQ.
    """
    required_fields = ['pbid', 'id']
    tally_by_loci = defaultdict(lambda: []) # PB.1 --> [(PB.1.1, zmw1), (PB.1.1, zmw2), (PB.1.2, zmw3)...]
    poor_zmws_not_in_rich = set()
    reader = DictReader(open(stat_filename),delimiter='\t')
    if any(x not in reader.fieldnames for x in required_fields):
        raise Exception("Expected fields `pbid`, `id` in {0} but only saw: {1}".format(\
            stat_filename, reader.fieldnames))

    for r in reader:
           m = rex_pbid.match(r['pbid'])
           if m is None:
               raise Exception("Expected PBID format PB.X.Y but saw {0}".format(r['pbid']))
           locus = m.group(1) # ex: PB.1
           m = rex_flnc.match(r['id'])
           if m is not None:
               zmw = m.group(1)
           else:
               zmw = r['id']
           if zmw in rich_zmws:
               tally_by_loci[locus].append((r['pbid'], zmw))
           else:
               poor_zmws_not_in_rich.add(zmw)

    return tally_by_loci, poor_zmws_not_in_rich

from collections import namedtuple
LocusInfo = namedtuple('LocusInfo', ['chrom', 'strand', 'regions', 'isoforms'])

def read_GFF(gff_filename, logf):
    """
    Read a GFF filename and get the gene regions

    :return: dict of (PB.X) --> LocusInfo
    """
    gff_info = {} # loci --> LocusInfo
    tmp = {} # loci PB.X --> list of GFF records for PB.X.Y

    for r in collapseGFFReader(gff_filename):
        m = rex_pbid.match(r.seqid)
        if m is None: raise Exception("Expected PBID format PB.X.Y but saw {0}".format(r.seqid))
        locus = m.group(1) # ex: PB.1
        if locus not in tmp:
            tmp[locus] = [r]
            gff_info[locus] = LocusInfo(chrom=r.chr, strand=r.strand, regions=None, isoforms=None)
        else:
            if gff_info[locus].chrom!=r.chr:
                logf.write("WARNING: Expected {0} to be on {1} but saw {2}. Could be minimap2 multi-mapping inconsistency for repetitive genes. Check later.\n".format(\
                    r.seqid, gff_info[locus].chrom, r.chr))
            tmp[locus].append(r)


    #pdb.set_trace()
    # now figure out the exonic regions for each gene PB.X
    for locus, records in tmp.items():
        c = ClusterTree(0, 0)
        for r in records:
            for e in r.ref_exons:
                c.insert(max(0,e.start-extra_bp_around_junctions), e.end+extra_bp_around_junctions, 1)

        regions = [(a,b) for (a,b,junk) in c.getregions()]
        regions[0] = (max(0,regions[0][0]-__padding_before_after__), regions[0][1])
        regions[-1] = (max(0,regions[-1][0]), regions[-1][1]+__padding_before_after__)
        gff_info[locus] = LocusInfo(chrom=gff_info[locus].chrom,
                                       strand=gff_info[locus].strand,
                                       regions=regions,
                                       isoforms=[r.seqid for r in records])

    return gff_info


def make_fake_genome(genome_d, gff_info, locus, output_prefix, output_name):

    chrom = gff_info[locus].chrom
    regions = gff_info[locus].regions

    with open(output_prefix+'.fasta', 'w') as f:
        f.write(">" + output_name + "\n")
        for s,e in regions:
            f.write(str(genome_d[chrom][s:e].seq))
        f.write("\n")
        f.close()

    # for mapping, write <0-based index on fake genome>, <ref chrom>, <0-based index on ref genome>
    with open(output_prefix+'.mapping.txt', 'w') as f:
        i = 0
        for s,e in regions:
            for j in range(s,e):
                f.write("{0},{1},{2}\n".format(i, chrom, j))
                i += 1

    with open(output_prefix+'.pbids.txt', 'w') as f:
        f.write("\n".join(gff_info[locus].isoforms)+'\n')

    print("Output written to {0}.fasta, {0}.mapping.txt, {0}.pbids.txt.".format(output_prefix), file=sys.stderr)


def select_loci_to_phase(args, genome_dict):

    logf = open('warning.logs', 'w')
    print("Reading FLNC file...", file=sys.stderr)
    flnc_fafq_d, rich_zmws = read_flnc_fafq(args.flnc_filename, 'fastq' if args.fq else 'fasta')

    print("Reading read_stat....", file=sys.stderr)
    tally_by_loci, poor_zmws_not_in_rich = read_read_stat(args.stat_filename, rich_zmws)

    print("Reading GFF file....", file=sys.stderr)
    gff_info = read_GFF(args.gff_filename, logf)

    # find all gene loci that has at least X FLNC coverage
    cand_loci = [k for k in tally_by_loci if len(tally_by_loci[k]) >= args.coverage]
    #pdb.set_trace()
    print("Total {0} loci read. {1} has >= {2} coverage.".format(\
        len(tally_by_loci), len(cand_loci), args.coverage), file=sys.stderr)
    for locus in cand_loci:
        if locus not in gff_info:
            logf.write("WARNING: {0} skipped because not in GFF info (probably filtered out later).\n".format(\
                locus))
            continue

        print("making", locus, file=sys.stderr)
        d2 = os.path.join("by_loci/{0}_size{1}".format(locus, len(tally_by_loci[locus])))
        os.makedirs(d2)

        chrom_len = dict((k, len(v)) for k,v in genome_dict.items())

        with open(os.path.join(d2, 'config'), 'w') as f:
            ref_start = max(0, gff_info[locus].regions[0][0])   # ref start must be >= 0
            # ref end must be at most chrom length
            ref_end = min(chrom_len[gff_info[locus].chrom], gff_info[locus].regions[-1][-1])

            #_chr, _strand, _start, _end = gff_info[locus]
            f.write("pbid={0}\n".format(locus))
            f.write("ref_chr={0}\n".format(gff_info[locus].chrom))
            f.write("ref_strand={0}\n".format(gff_info[locus].strand))
            f.write("ref_start={0}\n".format(ref_start))
            f.write("ref_end={0}\n".format(ref_end))

        make_fake_genome(genome_dict, gff_info, locus,
                         d2 + '/fake',
                         'fake_' + locus)


        # write ccs.fastq
        f1 = open(os.path.join(d2, 'ccs.fastq'), 'w')
        f2 = open(os.path.join(d2, 'ccs.fasta'), 'w')
        h = open(os.path.join(d2, 'fake.read_stat.txt'), 'w')
        h.write("id\tlength\tis_fl\tstat\tpbid\n")
        zmws_seen = set()
        for pbid, zmw in tally_by_loci[locus]:
            if zmw in zmws_seen: continue # duplicates can occassionally happen w split mapping
            zmws_seen.add(zmw)
            rec = flnc_fafq_d[zmw]
            if not args.fq: rec.letter_annotations['phred_quality'] = [60]*len(rec.seq) # manually addin QVs
            SeqIO.write(rec, f1, 'fastq')
            SeqIO.write(rec, f2, 'fasta')
            h.write("{0}\t{1}\tY\tunique\t{2}\n".format(zmw, len(rec.seq), pbid))
        f1.close()
        f2.close()
        h.close()
    logf.close()

def getargs():
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("genome_fasta", help="Reference genome fasta")
    parser.add_argument("flnc_filename", help="FLNC fasta/fastq file (if fastq, run with --fq)")
    parser.add_argument("gff_filename", help="GFF file of transcripts, IDs must be PB.X.Y")
    parser.add_argument("stat_filename", help="Tab-delimited read stat file linking FLNC to PB.X.Y")
    parser.add_argument("-c", "--coverage", type=int, default=40, help="Minimum FLNC coverage required (default: 40)")
    parser.add_argument("--fq", default=False, action="store_true", help="Use if flnc is fastq instead of fasta")

    return parser

if __name__ == "__main__":
    print("Reading genome...", file=sys.stderr)
    parser = getargs()

    args = parser.parse_args()

    if os.path.exists("by_loci"):
        print("Directory by_loci/ already exists. Delete before running!", file=sys.stderr)
        sys.exit(-1)

    if not os.path.exists(args.genome_fasta):
        print("Cannot find genome FASTA {0}. Abort!".format(args.genome_fasta), file=sys.stderr)
        sys.exit(-1)

    if not os.path.exists(args.flnc_filename):
        print("Cannot find FLNC file {0}. Abort!".format(args.flnc_filename), file=sys.stderr)
        sys.exit(-1)

    if not os.path.exists(args.gff_filename):
        print("Cannot find GFF file {0}. Abort!".format(args.gff_filename), file=sys.stderr)
        sys.exit(-1)

    if not os.path.exists(args.stat_filename):
        print("Cannot find Stat file {0}. Abort!".format(args.stat_filename), file=sys.stderr)
        sys.exit(-1)

    print("Reading genome fasta {0}....".format(args.genome_fasta), file=sys.stderr)
    genome_d = SeqIO.to_dict(SeqIO.parse(open(args.genome_fasta),'fasta'))

    select_loci_to_phase(args, genome_d)
