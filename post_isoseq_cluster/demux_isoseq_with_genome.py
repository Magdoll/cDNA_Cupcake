#!/usr/bin/env python
__author__ = 'etseng@pacb.com'

"""
Demultiplex IsoSeq1/IsoSeq2/IsoSeq3 job output (with genome mapping)

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
    (isoseq1 and isoseq2)
    id      length  is_fl   stat    pbid
    m54006_170729_232022/43123426/1712_71_CCS       1641    Y       unique  PB.3811.1
    m54006_170729_232022/26476826/32_1680_CCS       1648    Y       unique  PB.3811.1
    m54006_170729_232022/44958075/3189_67_CCS       3122    Y       unique  PB.3787.1
    m54006_170729_232022/9896496/3222_72_CCS        3150    Y       unique  PB.3787.1
    (isoseq3)
    id      length  is_fl   stat    pbid
    m54043_180729_193105/31392528/ccs       5306    Y       unique  PB.6318.2
    m54043_180729_193105/53346556/ccs       5233    Y       unique  PB.6318.2
    m54043_180729_193105/66716361/ccs       5030    Y       unique  PB.8713.2
    m54043_180729_193105/6488380/ccs        5001    Y       unique  PB.8713.2
    
Classify report format:
   (isoseq1 and 2)
    id,strand,fiveseen,polyAseen,threeseen,fiveend,polyAend,threeend,primer,chimera
    m54033_171031_152256/27919078/31_1250_CCS,+,1,1,1,31,1250,1280,3,0
    m54033_171031_152256/27919079/31_3840_CCS,+,1,1,1,31,3840,3869,2,0
    m54033_171031_152256/27919086/29_3644_CCS,+,1,1,1,29,3644,3674,2,0

    (isoseq3)
    id,strand,fivelen,threelen,polyAlen,insertlen,primer_index,primer
    m54020_170625_150952/69664969/ccs,-,31,39,57,2627,0--7,Clontech--bc7
    m54020_170625_150952/69664996/ccs,-,31,40,59,990,0--6,Clontech--bc6
    m54020_170625_150952/69664999/ccs,+,30,38,59,1724,0--1,Clontech--bc1
"""

import os, re, sys
from csv import DictReader, DictWriter
from collections import defaultdict, Counter
from Bio import SeqIO

hq1_id_rex = re.compile('i\d+_HQ_\S+\|(\S+)\/f\d+p\d+\/\d+')
hq2_id_rex = re.compile('HQ_\S+\|(\S+)\/f\d+p\d+\/\d+')
hq3_id_rex = re.compile('[\S_]+(transcript/\d+)')

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
    # location for mapped fastq in IsoSeq3
    mapped_fastq3 = os.path.join(os.path.abspath(src_dir), 'tasks', 'isocollapse.tasks.collapse_mapped_isoforms-0', 'out.fastq')
    mapped_gff3 = os.path.join(os.path.abspath(src_dir), 'tasks', 'isocollapse.tasks.collapse_mapped_isoforms-0', 'out.gff')
    # location for read_stat.txt in IsoSeq1 IsoSeq2
    read_stat = os.path.join(os.path.abspath(src_dir), 'tasks', 'pbtranscript.tasks.post_mapping_to_genome-0', 'output_mapped.no5merge.collapsed.read_stat.txt')
    read_stat2 = os.path.join(os.path.abspath(src_dir), 'tasks', 'pbtranscript2tools.tasks.post_mapping_to_genome-0', 'output_mapped.no5merge.collapsed.read_stat.txt')
    read_stat3 = os.path.join(os.path.abspath(src_dir), 'tasks', 'isocollapse.tasks.collapse_mapped_isoforms-0', 'out.no5merge.collapsed.read_stat.txt')
    # location for classify report in IsoSeq1 and 2
    primer_csv = os.path.join(os.path.abspath(src_dir), 'tasks', 'pbcoretools.tasks.gather_csv-1', 'file.csv')
    lima_report = os.path.join(os.path.abspath(src_dir), 'tasks', 'barcoding.tasks.lima-0', 'lima_output.lima.report')

    if os.path.exists(mapped_fastq):
        print >> sys.stderr, "Detecting IsoSeq1 task directories..."
        os.symlink(mapped_fastq, os.path.join(out_dir, 'mapped.fastq'))
        os.symlink(mapped_gff, os.path.join(out_dir, 'mapped.gff'))
        os.symlink(read_stat, os.path.join(out_dir, 'mapped.read_stat.txt'))
        os.symlink(primer_csv, os.path.join(out_dir, 'classify_report.csv'))
        isoseq_version = '1'
    elif os.path.exists(mapped_fastq2):
        print >> sys.stderr, "Detecting IsoSeq2 task directories..."
        os.symlink(mapped_fastq2, os.path.join(out_dir, 'mapped.fastq'))
        os.symlink(mapped_gff2, os.path.join(out_dir, 'mapped.gff'))
        os.symlink(read_stat2, os.path.join(out_dir, 'mapped.read_stat.txt'))
        os.symlink(primer_csv, os.path.join(out_dir, 'classify_report.csv'))
        isoseq_version = '2'
    elif os.path.exists(mapped_fastq3):
        print >> sys.stderr, "Detecting IsoSeq3 task directories..."
        os.symlink(mapped_fastq3, os.path.join(out_dir, 'mapped.fastq'))
        os.symlink(mapped_gff3, os.path.join(out_dir, 'mapped.gff'))
        os.symlink(read_stat3, os.path.join(out_dir, 'mapped.read_stat.txt'))
        make_classify_csv_from_lima_report(lima_report, os.path.join(out_dir, 'classify_report.csv'))
        isoseq_version = '3'
    else:
        print >> sys.stderr, "Cannot find HQ FASTQ in job directory! Does not look like Iso-Seq1, 2, or 3 jobs!"
        sys.exit(-1)
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
    :return: primer list, dict of FL/nFL id --> primer (in integer if isoseq1 or 2)
    """
    # for isoseq1 and 2, primer is an integer
    # for isoseq3, use primer_index field which is ex: 0--1, 0--2
    info = {}
    primer_list = set()
    for r in DictReader(open(classify_csv), delimiter=','):
        if r['primer']=='NA': continue # skip nFL
        if 'primer_index' in r: # isoseq3 version
            p = r['primer_index']
        else: # isoseq1 or 2
            p = r['primer']
        primer_list.add(p)
        info[r['id']] = p
    return primer_list, info

def make_classify_csv_from_lima_report(report_filename, output_filename):

    h = open(output_filename, 'w')
    fout = DictWriter(h, fieldnames=['id', 'primer_index', 'primer'], delimiter=',')
    fout.writeheader()
    for r in DictReader(open(report_filename), delimiter='\t'):
        assert r['IdxLowestNamed'].endswith('_5p') or r['IdxLowestNamed'].endswith('_3p')
        assert r['IdxHighestNamed'].endswith('_5p') or r['IdxHighestNamed'].endswith('_3p')
        if r['IdxLowestNamed'].endswith('_5p') and r['IdxHighestNamed'].endswith('_3p'):
            newrec = {'id': r['ZMW']+'/ccs',
                    'primer_index': "{0}--{1}".format(r['IdxLowest'], r['IdxHighest']),
                    'primer': "{0}--{1}".format(r['IdxLowestNamed'], r['IdxHighestNamed'])}
        fout.writerow(newrec)
    h.close()


def main(job_dir=None, mapped_fastq=None, read_stat=None, classify_csv=None, output_filename=sys.stdout, primer_names=None):
    if job_dir is not None:
        out_dir_ignore, mapped_fastq, read_stat, classify_csv, isoseq_version = link_files(job_dir)
        assert isoseq_version in ('1', '2', '3')
    else:
        assert os.path.exists(mapped_fastq)
        assert os.path.exists(read_stat)
        assert os.path.exists(classify_csv)

    # info: dict of hq_isoform --> primer --> FL count
    print >> sys.stderr, "Reading {0}....".format(classify_csv)
    primer_list, classify_csv = read_classify_csv(classify_csv)
    print >> sys.stderr, "Reading {0}....".format(read_stat)
    info = read_read_stat(read_stat, classify_csv)

    primer_list = list(primer_list)
    primer_list.sort()
    # if primer names are not given, just use as is...
    tmp_primer_names = dict((x,x) for x in primer_list)
    if primer_names is None:
        primer_names = tmp_primer_names
    else:
        for k,v in tmp_primer_names.iteritems():
            if k not in primer_names:
                primer_names[k] = v

    f = open(output_filename, 'w')
    f.write("id,{0}\n".format(",".join(primer_names.values())))
    print >> sys.stderr, "Reading {0}....".format(mapped_fastq)
    for r in SeqIO.parse(open(mapped_fastq), 'fastq'):
        pbid = r.id.split('|')[0]
        f.write(pbid)
        for p in primer_names:
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
            primer_names[index] = name
    else:
        primer_names = None

    main(args.job_dir, args.mapped_fastq, args.read_stat, args.classify_csv, args.output, primer_names)
