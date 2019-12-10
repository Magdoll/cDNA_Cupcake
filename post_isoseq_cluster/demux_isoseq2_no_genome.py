#!/usr/bin/env python
__author__ = 'etseng@pacb.com'

"""
Demultiplex IsoSeq1/IsoSeq2 job output (without genome mapping)

INPUT: HQ isoform and report.csv and file.csv (alternatively, job directory)
OUTPUT: CSV file containing associated FL count for each isoform

HQ isoform ID format:
  (isoseq1) i0_HQ_sample3f1db2|c44/f3p0/2324
  (isoseq2) HQ_sampleAZxhBguy|cb1063_c22/f3p0/6697
Cluster report format:
   (isoseq1)
    cluster_id,read_id,read_type
    i0_ICE_sample3f1db2|c23,m54033_171031_152256/26476965/30_5265_CCS,FL
    i0_ICE_sample3f1db2|c43,m54033_171031_152256/32441242/25_2283_CCS,FL
    i0_ICE_sample3f1db2|c43,m54033_171031_152256/50201527/30_2323_CCS,FL
    i0_ICE_sample3f1db2|c44,m54033_171031_152256/49545897/2374_68_CCS,FL
   (isoseq2)
    cluster_id,read_id,read_type
    cb10063_c6,m54006_170729_232022/56361426/29_9138_CCS,FL
    cb10407_c1,m54006_170729_232022/16712197/2080_72_CCS,FL
    cb10467_c49,m54006_170729_232022/48104064/0_2421_CCS,NonFL
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
    Locate HQ isoform, (cluster) report.csv, (classify) file.csv link to current directory
    """
    # location for HQ fastq in IsoSeq1
    hq_fastq = os.path.join(os.path.abspath(src_dir), 'tasks', 'pbtranscript.tasks.combine_cluster_bins-0', 'hq_isoforms.fastq')
    # location for HQ fastq in IsoSeq2
    hq_fastq2 = os.path.join(os.path.abspath(src_dir), 'tasks', 'pbtranscript2tools.tasks.collect_polish-0', 'all_arrowed_hq.fastq')
    # location for cluster report in IsoSeq1
    cluster_csv = os.path.join(os.path.abspath(src_dir), 'tasks', 'pbtranscript.tasks.combine_cluster_bins-0', 'cluster_report.csv')
    cluster_csv2 = os.path.join(os.path.abspath(src_dir), 'tasks', 'pbtranscript2tools.tasks.collect_polish-0', 'report.csv')
    # location for classify report in IsoSeq1 and 2
    primer_csv = os.path.join(os.path.abspath(src_dir), 'tasks', 'pbcoretools.tasks.gather_csv-1', 'file.csv')

    if os.path.exists(hq_fastq):
        print("Detecting IsoSeq1 task directories...", file=sys.stderr)
        os.symlink(hq_fastq, os.path.join(out_dir, 'hq_isoforms.fastq'))
        os.symlink(cluster_csv, os.path.join(out_dir, 'cluster_report.csv'))
        os.symlink(primer_csv, os.path.join(out_dir, 'classify_report.csv'))
        isoseq_version = '1'
    else:
        print("Detecting IsoSeq2 task directories...", file=sys.stderr)
        os.symlink(hq_fastq2, os.path.join(out_dir, 'hq_isoforms.fastq'))
        os.symlink(cluster_csv2, os.path.join(out_dir, 'cluster_report.csv'))
        os.symlink(primer_csv, os.path.join(out_dir, 'classify_report.csv'))
        isoseq_version = '2'
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
                cid = r['cluster_id'].split('|')[1]
            else:
                cid = r['cluster_id']
            info[cid][p] += 1
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


def main(job_dir=None, hq_fastq=None, cluster_csv=None, classify_csv=None, output_filename=sys.stdout):
    if job_dir is not None:
        out_dir_ignore, hq_fastq, cluster_csv, classify_csv, isoseq_version = link_files(job_dir)
        assert isoseq_version in ('1', '2')
    else:
        assert os.path.exists(hq_fastq)
        assert os.path.exists(cluster_csv)
        assert os.path.exists(classify_csv)

    # info: dict of hq_isoform --> primer --> FL count
    print("Reading {0}....".format(classify_csv), file=sys.stderr)
    max_primer, classify_csv = read_classify_csv(classify_csv)
    print("Reading {0}....".format(cluster_csv), file=sys.stderr)
    info = read_cluster_csv(cluster_csv, classify_csv, isoseq_version)

    f = open(output_filename, 'w')
    f.write("id,{0}\n".format(",".join("primer"+str(i) for i in range(max_primer+1))))
    print("Reading {0}....".format(hq_fastq), file=sys.stderr)
    for r in SeqIO.parse(open(hq_fastq), 'fastq'):
        if isoseq_version=='1':
            m = hq1_id_rex.match(r.id)
        else:
            m = hq2_id_rex.match(r.id)

        if m is None:
            print("Unexpected HQ isoform ID format: {0}! Abort.".format(r.id), file=sys.stderr)
            sys.exit(-1)
        cid = m.group(1)
        f.write(r.id)
        for p in range(max_primer+1):
            f.write(",{0}".format(info[cid][p]))
        f.write("\n")
    f.close()
    print("Count file written to {0}.".format(f.name), file=sys.stderr)


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser("")
    parser.add_argument("-j", "--job_dir", help="Job directory (if given, automatically finds required files)")
    parser.add_argument("--hq_fastq", help="HQ isoform fastq (overridden by --job_dir if given)")
    parser.add_argument("--cluster_csv", help="Cluster report CSV (overridden by --job_dir if given)")
    parser.add_argument("--classify_csv", help="Classify report CSV (overriden by --job_dir if given)")
    parser.add_argument("-o", "--output", help="Output count filename", required=True)

    args = parser.parse_args()


    main(args.job_dir, args.hq_fastq, args.cluster_csv, args.classify_csv, args.output)
