__author__ = "etseng@pacb.com"

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
from cupcake.io.SeqReaders import LazyFastqReader
from cupcake.io.GFF import collapseGFFReader

rex_flnc = re.compile('(m\d+_\d+_\d+\/\d+)\/ccs')
rex_pbid = re.compile('(PB.\d+).\d+')


extra_bp_around_junctions = 50 # get this much around junctions to be safe AND to not screw up GMAP who doesn't like microintrons....
__padding_before_after__ = 10 # get this much before and after the start


def read_flnc_fastq(flnc_filename):
    """
    Read FLNC fastq into a dict of zmw --> lazy file pointer
    """
    flnc_fastq_d = LazyFastqReader(flnc_filename)
    rich_zmws = set()

    for k in flnc_fastq_d.keys():
        m = rex_flnc.match(k)
        if m is None:
            raise Exception, "Expected FLNC id format is <movie>/<zmw>/ccs! Instead saw: {0}".format(k)
        zmw = m.group(1)
        flnc_fastq_d.d[zmw] = flnc_fastq_d.d[k]
        rich_zmws.add(zmw)

    return flnc_fastq_d, rich_zmws


def read_read_stat(stat_filename, rich_zmws):
    """
    Read .read_stat.txt file
    :return: tally_by_loci -- dict  of {locus -- list of (isoform, zmw)}
             poor_zmws_not_in_rich --- list of ZMWs that were in .read_stat but missing in FLNC FASTQ. ONly happens if used diff CCS to get FASTQ.
    """
    required_fields = ['pbid', 'id', 'is_fl', 'stat']
    tally_by_loci = defaultdict(lambda: []) # PB.1 --> [(PB.1.1, zmw1), (PB.1.1, zmw2), (PB.1.2, zmw3)...]
    poor_zmws_not_in_rich = set()
    reader = DictReader(open(stat_filename),delimiter='\t')
    if any(x not in reader.fieldnames for x in required_fields):
        raise Exception, "Expected fields `pbid`, `is_fl`, and `stat` in {0} but only saw: {1}".format(\
            stat_filename, reader.fieldnames)

    for r in reader:
       if r['is_fl']=='Y' and r['stat']=='unique':
           m = rex_pbid.match(r['pbid'])
           if m is None:
               raise Exception, "Expected PBID format PB.X.Y but saw {0}".format(r['pbid'])
           locus = m.group(1) # ex: PB.1
           m = rex_flnc.match(r['id'])
           if m is None:
               raise Exception, "Expected FLNC id format is <movie>/<zmw>/ccs! Instead saw: {0}".format(r['id'])
           zmw = m.group(1)
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
        if m is None: raise Exception, "Expected PBID format PB.X.Y but saw {0}".format(r.seqid)
        locus = m.group(1) # ex: PB.1
        if locus not in tmp:
            tmp[locus] = [r]
            gff_info[locus] = LocusInfo(chrom=r.chr, strand=r.strand, regions=None, isoforms=None)
        else:
            if gff_info[locus].chrom!=r.chr:
                logf.write("WARNING: Expected {0} to be on {1} but saw {2}. Could be minimap2 multi-mapping inconsistency for repetitive genes. Check later.\n".format(\
                    r.seqid, gff_info[locus].chrom, r.chr))
            tmp[locus].append(r)


    # now figure out the exonic regions for each gene PB.X
    for locus, records in tmp.iteritems():
        c = ClusterTree(0, 0)
        for r in records:
            for e in r.ref_exons:
                c.insert(e.start-extra_bp_around_junctions, e.end+extra_bp_around_junctions, 1)

        regions = [(a,b) for (a,b,junk) in c.getregions()]
        regions[0] = (regions[0][0]-__padding_before_after__, regions[0][1])
        regions[-1] = (regions[-1][0], regions[-1][1]+__padding_before_after__)
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
            for j in xrange(s,e):
                f.write("{0},{1},{2}\n".format(i, chrom, j))
                i += 1

    with open(output_prefix+'.pbids.txt', 'w') as f:
        f.write("\n".join(gff_info[locus].isoforms)+'\n')

    print >> sys.stderr, "Output written to {0}.fasta, {0}.mapping.txt, {0}.pbids.txt.".format(output_prefix)


def select_loci_to_phase(args, genome_dict):

    logf = open('warning.logs', 'w')
    print >> sys.stderr, "Reading FLNC file..."
    flnc_fastq_d, rich_zmws = read_flnc_fastq(args.flnc_filename)

    print >> sys.stderr, "Reading read_stat...."
    tally_by_loci, poor_zmws_not_in_rich = read_read_stat(args.stat_filename, rich_zmws)

    print >> sys.stderr, "Reading GFF file...."
    gff_info = read_GFF(args.gff_filename, logf)

    # find all gene loci that has at least X FLNC coverage
    cand_loci = filter(lambda k: len(tally_by_loci[k]) >= args.coverage, tally_by_loci)
    print >> sys.stderr, "Total {0} loci read. {1} has >= {2} coverage.".format(\
        len(tally_by_loci), len(cand_loci), args.coverage)
    for locus in cand_loci:
        if locus not in gff_info:
            logf.write("WARNING: {0} skipped because not in GFF info (probably filtered out later).\n".format(\
                locus))
            continue

        print >> sys.stderr, "making", locus
        d2 = os.path.join("by_loci/{0}_size{1}".format(locus, len(tally_by_loci[locus])))
        os.makedirs(d2)
        with open(os.path.join(d2, 'config'), 'w') as f:
            #_chr, _strand, _start, _end = gff_info[locus]
            f.write("pbid={0}\n".format(locus))
            f.write("ref_chr={0}\n".format(gff_info[locus].chrom))
            f.write("ref_strand={0}\n".format(gff_info[locus].strand))
            f.write("ref_start={0}\n".format(gff_info[locus].regions[0][0]))
            f.write("ref_end={0}\n".format(gff_info[locus].regions[-1][-1]))

        make_fake_genome(genome_dict, gff_info, locus,
                         d2 + '/fake',
                         'fake_' + locus)


        # write ccs.fastq
        f = open(os.path.join(d2, 'ccs.fastq'), 'w')
        h = open(os.path.join(d2, 'fake.read_stat.txt'), 'w')
        h.write("id\tlength\tis_fl\tstat\tpbid\n")
        for pbid, zmw in tally_by_loci[locus]:
            rec = flnc_fastq_d[zmw]
            SeqIO.write(rec, f, 'fastq')
            h.write("{0}\t{1}\tY\tunique\t{2}\n".format(zmw, len(rec.seq), pbid))
        f.close()
        h.close()
    logf.close()

def getargs():
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("genome_fasta")
    parser.add_argument("flnc_filename")
    parser.add_argument("gff_filename")
    parser.add_argument("stat_filename")
    parser.add_argument("-c", "--coverage", type=int, default=40, help="Minimum FLNC coverage required (default: 40)")

    return parser

if __name__ == "__main__":
    print >> sys.stderr, "Reading genome..."
    genome_d = SeqIO.to_dict(SeqIO.parse(open('B73_RefV4.fa'),'fasta'))
