#!/usr/bin/env python
__author__ = 'etseng@pacb.com'
import pdb
"""
Demultiplex IsoSeq (SMRT Link 8.0) job output (with genome mapping)
"""

import os, re, sys
from csv import DictReader, DictWriter
from collections import defaultdict, Counter
from Bio import SeqIO

mapped_id_rex = re.compile('(PB\S*.\d+[.\d+]?)')

def type_fafq(fafq):
    x = fafq.upper()
    if x.endswith('.FA') or x.endswith('.FASTA'): return 'fasta'
    elif x.endswith('.FQ') or x.endswith('.FASTQ'): return 'fastq'
    else: raise Exception("Mapped fasta/fastq filename must end with .fasta or .fastq! Saw {0} instead, abort!".format(fafq))

def link_files(src_dir, out_dir='./'):
    """
    :param src_dir: job directory
    Locate mapped.fastq, read-stat, classify report link to current directory
    """
    # location for mapped fastq in IsoSeq3
    mapped_fastq = os.path.join(os.path.abspath(src_dir), 'outputs', 'collapse_isoforms.fastq')  # for <SL8
    mapped_fasta = os.path.join(os.path.abspath(src_dir), 'outputs', 'collapse_isoforms.fasta')  # SL8+ only fasta
    mapped_gff = os.path.join(os.path.abspath(src_dir), 'outputs', 'collapse_isoforms.gff')
    read_stat = os.path.join(os.path.abspath(src_dir), 'outputs', 'collapse_isoforms.read_stat.txt')
    primer_csv = os.path.join(os.path.abspath(src_dir), 'outputs', 'flnc.report.csv')

    if os.path.exists(mapped_fastq):
        print("Detecting IsoSeq task directories...", file=sys.stderr)
        return out_dir, mapped_fastq, read_stat, primer_csv
    elif os.path.exists(mapped_fasta):
        print("Detecting IsoSeq task directories...", file=sys.stderr)
        return out_dir, mapped_fasta, read_stat, primer_csv
    else:
        print("Cannot find expected files (ex: collapse_isoforms.fastq) in job directory! Does not look like a Iso-Seq job!", file=sys.stderr)
        sys.exit(-1)


def read_read_stat(read_stat, classify_info):
    """
    :return: dict of pbid --> (int) primer --> FL count
    """
    info = defaultdict(lambda: Counter())
    for r in DictReader(open(read_stat), delimiter='\t'):
        p = classify_info[r['id']]
        info[r['pbid']][p] += 1
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


def main(job_dir=None, mapped_fafq=None, read_stat=None, classify_csv=None, output_filename=sys.stdout, primer_names=None):
    if job_dir is not None:
        out_dir_ignore, mapped_fafq, read_stat, classify_csv = link_files(job_dir)
    else:
        assert os.path.exists(mapped_fafq)
        assert os.path.exists(read_stat)
        assert os.path.exists(classify_csv)

    # info: dict of hq_isoform --> primer --> FL count
    print("Reading {0}....".format(classify_csv), file=sys.stderr)
    primer_list, classify_info = read_classify_csv(classify_csv)
    print("Reading {0}....".format(read_stat), file=sys.stderr)
    info = read_read_stat(read_stat, classify_info)

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
    f.write("id,{0}\n".format(",".join(list(primer_names.values()))))
    print("Reading {0}....".format(mapped_fafq), file=sys.stderr)
    for r in SeqIO.parse(open(mapped_fafq), type_fafq(mapped_fafq)):
        m = mapped_id_rex.match(r.id)  # expected ID: PB.X.Y|xxxx.....
        if m is None:
            raise Exception("Expected ID format PB.X.Y but found {0}!".format(r.id))
        pbid = m.group(1)
        if pbid not in info:
            print("WARNING: Could not find {0} in {1}. No demux output for this sequence!".format(pbid, read_stat))
        else:
            #pdb.set_trace()
            f.write(pbid)
            for p in primer_names:
                f.write(",{0}".format(info[pbid][p]))
            f.write("\n")
    f.close()
    print("Count file written to {0}.".format(f.name), file=sys.stderr)


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser("")
    parser.add_argument("-j", "--job_dir", help="Job directory (if given, automatically finds required files)")
    parser.add_argument("--mapped_fafq", help="mapped fasta/fastq (overridden by --job_dir if given)")
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

    main(args.job_dir, args.mapped_fafq, args.read_stat, args.classify_csv, args.output, primer_names)
