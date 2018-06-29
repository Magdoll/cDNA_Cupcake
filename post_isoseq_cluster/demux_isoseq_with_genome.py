#!/usr/bin/env python
__author__ = 'etseng@pacb.com'

"""
Demultiplex IsoSeq1/IsoSeq2 job output (with genome mapping)

=================
WITHOUT GENOME
=================
INPUT: hq fastq, cluster report, classify report (alternatively, job directory)
OUTPUT: CSV file containing associated FL count for each isoform:primer

=================
WITH GENOME
=================
INPUT: out_mapped.fastq, read_stat.txt, classify report (alternatively, job directory)
OUTPUT: CSV file containing associated FL count for each isoform:primer

Read stat format:
    id      length  is_fl   stat    pbid
    m54006_170729_232022/43123426/1712_71_CCS       1641    Y       unique  PB.3811.1
    m54006_170729_232022/26476826/32_1680_CCS       1648    Y       unique  PB.3811.1
    m54006_170729_232022/44958075/3189_67_CCS       3122    Y       unique  PB.3787.1
    m54006_170729_232022/9896496/3222_72_CCS        3150    Y       unique  PB.3787.1

Classify report format:
   (isoseq1 and 2)
    id,strand,fiveseen,polyAseen,threeseen,fiveend,polyAend,threeend,primer,chimera
    m54033_171031_152256/27919078/31_1250_CCS,+,1,1,1,31,1250,1280,3,0
    m54033_171031_152256/27919079/31_3840_CCS,+,1,1,1,31,3840,3869,2,0
    m54033_171031_152256/27919086/29_3644_CCS,+,1,1,1,29,3644,3674,2,0
"""

import os, re, sys
from csv import DictReader
from collections import defaultdict, Counter
from Bio import SeqIO

hq1_id_rex = re.compile('i\d+_HQ_\S+\|(\S+)\/f\d+p\d+\/\d+')
hq2_id_rex = re.compile('HQ_\S+\|(\S+)\/f\d+p\d+\/\d+')

def link_files(src_dir, out_dir='./'):
    """
    :param src_dir: job directory
    Locate mapped.fastq, read-stat, classify report link to current directory
    """
    # location for mapped fastq in IsoSeq1
    mapped_fastq = os.path.join(os.path.abspath(src_dir), 'tasks', 'pbtranscript.tasks.post_mapping_to_genome-0', 'output_mapped.fastq')
    mapped_gff = os.path.join(os.path.abspath(src_dir), 'tasks', 'pbtranscript.tasks.post_mapping_to_genome-0', 'output_mapped.gff')
    # location for mapped fastq in IsoSeq2
    mapped_fastq2 = os.path.join(os.path.abspath(src_dir), 'tasks', 'pbtranscript2tools.tasks.post_mapping_to_genome-0', 'output_mapped.fastq')
    mapped_gff2 = os.path.join(os.path.abspath(src_dir), 'tasks', 'pbtranscript2tools.tasks.post_mapping_to_genome-0', 'output_mapped.gff')
    # location for read_stat.txt in IsoSeq1 IsoSeq2
    read_stat = os.path.join(os.path.abspath(src_dir), 'tasks', 'pbtranscript.tasks.post_mapping_to_genome-0', 'output_mapped.no5merge.collapsed.read_stat.txt')
    read_stat2 = os.path.join(os.path.abspath(src_dir), 'tasks', 'pbtranscript2tools.tasks.post_mapping_to_genome-0', 'output_mapped.no5merge.collapsed.read_stat.txt')
    # location for classify report in IsoSeq1 and 2
    primer_csv = os.path.join(os.path.abspath(src_dir), 'tasks', 'pbcoretools.tasks.gather_csv-1', 'file.csv')

    if os.path.exists(mapped_fastq):
        print >> sys.stderr, "Detecting IsoSeq1 task directories..."
        os.symlink(mapped_fastq, os.path.join(out_dir, 'mapped.fastq'))
        os.symlink(mapped_gff, os.path.join(out_dir, 'mapped.gff'))
        os.symlink(read_stat, os.path.join(out_dir, 'mapped.read_stat.txt'))
        os.symlink(primer_csv, os.path.join(out_dir, 'classify_report.csv'))
        isoseq_version = '1'
    else:
        print >> sys.stderr, "Detecting IsoSeq2 task directories..."
        os.symlink(mapped_fastq2, os.path.join(out_dir, 'mapped.fastq'))
        os.symlink(mapped_gff2, os.path.join(out_dir, 'mapped.gff'))
        os.symlink(read_stat2, os.path.join(out_dir, 'mapped.read_stat.txt'))
        os.symlink(primer_csv, os.path.join(out_dir, 'classify_report.csv'))
        isoseq_version = '2'
    return out_dir, 'mapped.fastq', 'mapped.read_stat.txt', 'classify_report.csv', isoseq_version

def read_read_stat(read_stat, classify_info):
    """
    :return: dict of pbid --> (int) primer --> FL count
    """
    info = defaultdict(lambda: Counter())
    for r in DictReader(open(read_stat), delimiter='\t'):
        if r['is_fl'] == 'Y':
            p = classify_info[r['id']]
            info[r['pbid']][p] += 1
    return dict(info)

def read_classify_csv(classify_csv):
    """
    :param classify_csv: classify report csv
    :return: primer range, dict of FL/nFL id --> primer (in integer)
    """
    info = {}
    max_p = 0
    for r in DictReader(open(classify_csv), delimiter=','):
        if r['primer']=='NA': continue # skip nFL
        p = int(r['primer'])
        max_p = max(max_p, p)
        info[r['id']] = p
    return max_p, info


def main(job_dir=None, mapped_fastq=None, read_stat=None, classify_csv=None, output_filename=sys.stdout, primer_names=None):
    if job_dir is not None:
        out_dir_ignore, mapped_fastq, read_stat, classify_csv, isoseq_version = link_files(job_dir)
        assert isoseq_version in ('1', '2')
    else:
        assert os.path.exists(mapped_fastq)
        assert os.path.exists(read_stat)
        assert os.path.exists(classify_csv)

    # info: dict of hq_isoform --> primer --> FL count
    print >> sys.stderr, "Reading {0}....".format(classify_csv)
    max_primer, classify_csv = read_classify_csv(classify_csv)
    print >> sys.stderr, "Reading {0}....".format(read_stat)
    info = read_read_stat(read_stat, classify_csv)

    # if primer names are not given, just use primer0, primer1, primer2, primer3....
    tmp_primer_names = dict((i,"primer"+str(i)) for i in xrange(max_primer+1))
    if primer_names is None:
        primer_names = tmp_primer_names
    else:
        for k,v in tmp_primer_names.iteritems():
            if k not in primer_names:
                primer_names[k] = v

    f = open(output_filename, 'w')
    f.write("id,{0}\n".format(",".join(primer_names[i] for i in xrange(max_primer+1))))
    print >> sys.stderr, "Reading {0}....".format(mapped_fastq)
    for r in SeqIO.parse(open(mapped_fastq), 'fastq'):
        pbid = r.id.split('|')[0]
        f.write(pbid)
        if pbid not in info:
            print >> sys.stderr, "WARNING: {0} not seen in .read_stat.txt!".format(pbid)
            for p in xrange(max_primer+1): f.write(",0")
        else:
            for p in xrange(max_primer+1):
                f.write(",{0}".format(info[pbid][p]))
        f.write("\n")
    f.close()
    print >> sys.stderr, "Count file written to {0}.".format(f.name)


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser("")
    parser.add_argument("-j", "--job_dir", help="Job directory (if given, automatically finds required files)")
    parser.add_argument("--mapped_fastq", help="mapped fastq (overridden by --job_dir if given)")
    parser.add_argument("--read_stat", help="read_stat txt (overridden by --job_dir if given)")
    parser.add_argument("--classify_csv", help="Classify report CSV (overriden by --job_dir if given)")
    parser.add_argument("--primer_names", default=None, help="Text file showing primer sample names (default: None)")
    parser.add_argument("-o", "--output", help="Output count filename", required=True)

    args = parser.parse_args()

    if args.primer_names is not None:
        primer_names = {}
        for line in open(args.primer_names):
            index, name = line.strip().split()
            primer_names[int(index)] = name
    else:
        primer_names = None

    main(args.job_dir, args.mapped_fastq, args.read_stat, args.classify_csv, args.output, primer_names)
