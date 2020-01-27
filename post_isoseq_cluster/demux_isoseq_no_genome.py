#!/usr/bin/env python
__author__ = 'etseng@pacb.com'

"""
Demultiplex IsoSeq (SMRT Link 8.0) job output (without genome mapping)
"""

import os, re, sys
from csv import DictReader, DictWriter
from collections import defaultdict, Counter
from Bio import SeqIO

#hq1_id_rex = re.compile('(i\d+_HQ_\S+\|\S+)\/f\d+p\d+\/\d+')
#hq2_id_rex = re.compile('HQ_\S+\|(\S+)\/f\d+p\d+\/\d+')
hq3_id_rex = re.compile('[\S]+(transcript/\d+)')  # ex: UnnamedSample_HQ_transcript/0

def type_fafq(fafq):
    x = fafq.upper()
    if x.endswith('.FA') or x.endswith('.FASTA'): return 'fasta'
    elif x.endswith('.FQ') or x.endswith('.FASTQ'): return 'fastq'
    else: raise Exception("HQ fasta/fastq filename must end with .fasta or .fastq! Saw {0} instead, abort!".format(fafq))

def link_files(src_dir, out_dir='./'):
    """
    :param src_dir: job directory
    Locate HQ isoform, (cluster) report.csv, (classify) file.csv link to current directory
    """
    # location for mapped fastq in IsoSeq3
    hq_fasta = os.path.join(os.path.abspath(src_dir), 'outputs', 'hq_transcripts.fasta')
    hq_fastq = os.path.join(os.path.abspath(src_dir), 'outputs', 'hq_transcripts.fastq')
    primer_csv = os.path.join(os.path.abspath(src_dir), 'outputs', 'flnc.report.csv')
    cluster_csv = os.path.join(os.path.abspath(src_dir), 'outputs', 'polished.cluster_report.csv')

    if os.path.exists(hq_fasta) or os.path.exists(hq_fastq):
        print("Detecting IsoSeq directories...", file=sys.stderr)
    else:
        print("Cannot find hq_transcripts.fasta/fastq in job directory! Does not look like SMRTLink 8 Iso-Seq job!", file=sys.stderr)
        sys.exit(-1)
    return out_dir, hq_fastq if os.path.exists(hq_fastq) else hq_fasta, cluster_csv, primer_csv

def read_cluster_csv(cluster_csv, classify_info):
    """
    :param report_csv: cluster_report.csv
    :return: dict of cluster_id --> integer primer --> FL count
    """
    info = defaultdict(lambda: Counter())
    for r in DictReader(open(cluster_csv), delimiter=','):
        assert r['read_type'] == 'FL'  # always FL for isoseq3
        p = classify_info[r['read_id']]
        cid = r['cluster_id']
        info[cid][p] += 1
    return dict(info)


def read_classify_csv(classify_csv):
    """
    :param classify_csv: classify report csv
    :return: primer list, dict of FL id --> primer
    """
    info = {}
    primer_list = set()
    for r in DictReader(open(classify_csv), delimiter=','):
        p = r['primer']
        primer_list.add(p)
        if r['id'] in info:
            raise Exception("{0} listed more than once in {1}!".format(r['id'], classify_csv))
        info[r['id']] = p
    return primer_list, info

def main(job_dir=None, hq_fafq=None, cluster_csv=None, classify_csv=None, output_filename=sys.stdout, primer_names=None):
    if job_dir is not None:
        out_dir_ignore, hq_fafq, cluster_csv, classify_csv = link_files(job_dir)
    else:
        assert os.path.exists(hq_fafq)
        assert os.path.exists(cluster_csv)
        assert os.path.exists(classify_csv)

    # info: dict of hq_isoform --> primer --> FL count
    print("Reading {0}....".format(classify_csv), file=sys.stderr)
    primer_list, classify_csv = read_classify_csv(classify_csv)
    print("Reading {0}....".format(cluster_csv), file=sys.stderr)
    info = read_cluster_csv(cluster_csv, classify_csv)

    primer_list = list(primer_list)
    primer_list.sort()
    # if primer names are not given, just use as is...
    tmp_primer_names = dict((x,x) for x in primer_list)
    if primer_names is None:
        primer_names = tmp_primer_names
    else:
        for k,v in tmp_primer_names.items():
            if k not in primer_names:
                primer_names[k] = v

    f = open(output_filename, 'w')
    f.write("id,{0}\n".format(",".join(list(primer_names.keys()))))
    print("Reading {0}....".format(hq_fafq), file=sys.stderr)
    for r in SeqIO.parse(open(hq_fafq), type_fafq(hq_fafq)):
        f.write(r.id)
        m = hq3_id_rex.match(r.id)
        cid = m.group(1)
        for p in primer_names:
            f.write(",{0}".format(info[cid][p]))
        f.write("\n")
    f.close()
    print("Count file written to {0}.".format(f.name), file=sys.stderr)


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser("")
    parser.add_argument("-j", "--job_dir", help="Job directory (if given, automatically finds required files)")
    parser.add_argument("--hq_fafq", help="HQ isoform fasta/fastq (overridden by --job_dir if given)")
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

    main(args.job_dir, args.hq_fafq, args.cluster_csv, args.classify_csv, args.output, primer_names)
