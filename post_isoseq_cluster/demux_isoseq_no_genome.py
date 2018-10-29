#!/usr/bin/env python
__author__ = 'etseng@pacb.com'

"""
Demultiplex IsoSeq1/IsoSeq2/IsoSeq3 job output (without genome mapping)

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


HQ isoform ID format:
  (isoseq1) i0_HQ_sample3f1db2|c44/f3p0/2324
  (isoseq2) HQ_sampleAZxhBguy|cb1063_c22/f3p0/6697
  (isoseq3) <biosample>_HQ_transcript/0
Cluster report format:
   (isoseq1)
    cluster_id,read_id,read_type
    i0_ICE_sample3f1db2|c23,m54033_171031_152256/26476965/30_5265_CCS,FL
    i0_ICE_sample3f1db2|c43,m54033_171031_152256/32441242/25_2283_CCS,FL
   (isoseq2)
    cluster_id,read_id,read_type
    cb10063_c6,m54006_170729_232022/56361426/29_9138_CCS,FL
    cb10467_c49,m54006_170729_232022/48104064/0_2421_CCS,NonFL
   (isoseq3)
    cluster_id,read_id,read_type
    transcript/16170,m54043_180729_193105/20644507/ccs,FL
    transcript/16170,m54043_180729_193105/19595761/ccs,FL
Classify report format:
   (isoseq1 and 2)
    id,strand,fiveseen,polyAseen,threeseen,fiveend,polyAend,threeend,primer,chimera
    m54033_171031_152256/27919078/31_1250_CCS,+,1,1,1,31,1250,1280,3,0
    m54033_171031_152256/27919079/31_3840_CCS,+,1,1,1,31,3840,3869,2,0
    m54033_171031_152256/27919086/29_3644_CCS,+,1,1,1,29,3644,3674,2,0
"""

import os, re, sys
from csv import DictReader, DictWriter
from collections import defaultdict, Counter
from Bio import SeqIO

hq1_id_rex = re.compile('(i\d+_HQ_\S+\|\S+)\/f\d+p\d+\/\d+')
hq2_id_rex = re.compile('HQ_\S+\|(\S+)\/f\d+p\d+\/\d+')
hq3_id_rex = re.compile('[\S_]+(transcript/\d+)')

def link_files(src_dir, out_dir='./'):
    """
    :param src_dir: job directory
    Locate HQ isoform, (cluster) report.csv, (classify) file.csv link to current directory
    """
    # location for HQ fastq in IsoSeq1
    hq_fastq = os.path.join(os.path.abspath(src_dir), 'tasks', 'pbtranscript.tasks.combine_cluster_bins-0', 'hq_isoforms.fastq')
    # location for HQ fastq in IsoSeq2
    hq_fastq2 = os.path.join(os.path.abspath(src_dir), 'tasks', 'pbtranscript2tools.tasks.collect_polish-0', 'all_arrowed_hq.fastq')
    # location for HQ fastq in IsoSeq3
    hq_fastq3 = os.path.join(os.path.abspath(src_dir), 'tasks', 'pbcoretools.tasks.bam2fastq_transcripts-0', 'hq_transcripts.fastq')
    # location for cluster report
    cluster_csv = os.path.join(os.path.abspath(src_dir), 'tasks', 'pbtranscript.tasks.combine_cluster_bins-0', 'cluster_report.csv')
    cluster_csv2 = os.path.join(os.path.abspath(src_dir), 'tasks', 'pbtranscript2tools.tasks.collect_polish-0', 'report.csv')
    cluster_csv3 = os.path.join(os.path.abspath(src_dir), 'tasks', 'pbcoretools.tasks.gather_csv-1', 'file.csv')
    # location for classify report in IsoSeq1 and 2
    primer_csv = os.path.join(os.path.abspath(src_dir), 'tasks', 'pbcoretools.tasks.gather_csv-1', 'file.csv')
    lima_report = os.path.join(os.path.abspath(src_dir), 'tasks', 'barcoding.tasks.lima-0', 'lima_output.lima.report')

    if os.path.exists(hq_fastq):
        print >> sys.stderr, "Detecting IsoSeq1 task directories..."
        os.symlink(hq_fastq, os.path.join(out_dir, 'hq_isoforms.fastq'))
        os.symlink(cluster_csv, os.path.join(out_dir, 'cluster_report.csv'))
        os.symlink(primer_csv, os.path.join(out_dir, 'classify_report.csv'))
        isoseq_version = '1'
    elif os.path.exists(hq_fastq2):
        print >> sys.stderr, "Detecting IsoSeq2 task directories..."
        os.symlink(hq_fastq2, os.path.join(out_dir, 'hq_isoforms.fastq'))
        os.symlink(cluster_csv2, os.path.join(out_dir, 'cluster_report.csv'))
        os.symlink(primer_csv, os.path.join(out_dir, 'classify_report.csv'))
        isoseq_version = '2'
    elif os.path.exists(hq_fastq3):
        print >> sys.stderr, "Detecting IsoSeq3 directories..."
        os.symlink(hq_fastq3, os.path.join(out_dir, 'hq_isoforms.fastq'))
        os.symlink(cluster_csv3, os.path.join(out_dir, 'cluster_report.csv'))
        print >> sys.stderr, "Making classify_report.csv because not yet in job directories..."
        make_classify_csv_from_lima_report(lima_report, os.path.join(out_dir, 'classify_report.csv'))
        isoseq_version = '3'
    else:
        print >> sys.stderr, "Cannot find HQ FASTQ in job directory! Does not look like Iso-Seq1, 2, or 3 jobs!"
        sys.exit(-1)
    return out_dir, 'hq_isoforms.fastq', 'cluster_report.csv', 'classify_report.csv', isoseq_version

def read_cluster_csv(cluster_csv, classify_info, isoseq_version):
    """
    :param report_csv: cluster_report.csv
    :return: dict of cluster_id --> integer primer --> FL count
    """
    info = defaultdict(lambda: Counter())
    for r in DictReader(open(cluster_csv), delimiter=','):
        if r['read_type'] == 'FL':
            p = classify_info[r['read_id']]
            if isoseq_version=='1':
                cid = r['cluster_id'].replace('_ICE_', '_HQ_')
            else:
                cid = r['cluster_id']
            info[cid][p] += 1
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
        #if 'primer_index' in r: # isoseq3 version
        #    p = r['primer_index']
        #else: # isoseq1 or 2
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

def main(job_dir=None, hq_fastq=None, cluster_csv=None, classify_csv=None, output_filename=sys.stdout, primer_names=None):
    if job_dir is not None:
        out_dir_ignore, hq_fastq, cluster_csv, classify_csv, isoseq_version = link_files(job_dir)
        assert isoseq_version in ('1', '2', '3')
    else:
        assert os.path.exists(hq_fastq)
        assert os.path.exists(cluster_csv)
        assert os.path.exists(classify_csv)
        # must figure out the isoseq_version automatically
        #HQ isoform ID format:
        #(isoseq1) i0_HQ_sample3f1db2|c44/f3p0/2324
        #(isoseq2) HQ_sampleAZxhBguy|cb1063_c22/f3p0/6697
        #(isoseq3) transcript/0
        r = SeqIO.parse(open(hq_fastq), 'fastq').next()
        m = hq1_id_rex.match(r.id)
        if m is not None:
            isoseq_version = '1'
            print >> sys.stderr, "IsoSeq1 ID format detected."
        else:
            m = hq2_id_rex.match(r.id)
            if m is not None:
                isoseq_version = '2'
                print >> sys.stderr, "IsoSeq2 ID format detected."
            else:
                m = hq3_id_rex.match(r.id)
                if m is not None:
                    isoseq_version = '3'
                    print >> sys.stderr, "IsoSeq3 ID format detected."
                else:
                    raise Exception, "Unrecognized HQ sequence ID format for {0}!".format(r.id)

    # info: dict of hq_isoform --> primer --> FL count
    print >> sys.stderr, "Reading {0}....".format(classify_csv)
    primer_list, classify_csv = read_classify_csv(classify_csv)
    print >> sys.stderr, "Reading {0}....".format(cluster_csv)
    info = read_cluster_csv(cluster_csv, classify_csv, isoseq_version)

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
    f.write("id,{0}\n".format(",".join(primer_names.keys())))
    print >> sys.stderr, "Reading {0}....".format(hq_fastq)
    for r in SeqIO.parse(open(hq_fastq), 'fastq'):
        if isoseq_version=='1':
            m = hq1_id_rex.match(r.id)
        elif isoseq_version=='2':
            m = hq2_id_rex.match(r.id)
        else: # isoseq_version=3
            m = hq3_id_rex.match(r.id)

        if m is None:
            print >> sys.stderr, "Unexpected HQ isoform ID format: {0}! Abort.".format(r.id)
            sys.exit(-1)
        cid = m.group(1)
        f.write(r.id)
        for p in primer_names:
            f.write(",{0}".format(info[cid][p]))
        f.write("\n")
    f.close()
    print >> sys.stderr, "Count file written to {0}.".format(f.name)


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser("")
    parser.add_argument("-j", "--job_dir", help="Job directory (if given, automatically finds required files)")
    parser.add_argument("--hq_fastq", help="HQ isoform fastq (overridden by --job_dir if given)")
    parser.add_argument("--cluster_csv", help="Cluster report CSV (overridden by --job_dir if given)")
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

    main(args.job_dir, args.hq_fastq, args.cluster_csv, args.classify_csv, args.output, primer_names)
