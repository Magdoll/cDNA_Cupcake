__author__ = 'etseng@pacb.com'
__version__ = '1.2'

import os, sys, subprocess
from csv import DictReader
from collections import defaultdict
"""
pbtranscript.tasks.combine_cluster_bins-0
"""

GMAP_BIN = "/home/UNIXHOME/etseng/bin/gmap"
GMAP_DB = "/home/UNIXHOME/etseng/share/gmap_db_new/"
GMAP_CPUS = 12
SIRV_DIR = "/home/UNIXHOME/etseng/projects2016/Lexogen_SIRV/ground_truth/"

def link_files(src_dir, out_dir):
    os.makedirs(out_dir)
    hq_fastq = os.path.join(os.path.abspath(src_dir), 'tasks', 'pbtranscript.tasks.combine_cluster_bins-0', 'hq_isoforms.fastq')
    cluster_csv = os.path.join(os.path.abspath(src_dir), 'tasks', 'pbtranscript.tasks.combine_cluster_bins-0', 'cluster_report.csv')

    os.symlink(hq_fastq, os.path.join(out_dir, 'hq_isoforms.fastq'))
    os.symlink(cluster_csv, os.path.join(out_dir, 'cluster_report.csv'))

    return out_dir, 'hq_isoforms.fastq', 'cluster_report.csv'

def make_abundance_from_Sequel_cluster_csv(cluster_csv, collapse_prefix):
    """
    cluster_id,read_id,read_type
    i0_ICE_samplee5686f|c1,m54011_160718_221804/53150541/967_60_CCS,FL
    """
    fl_ass = defaultdict(lambda: set()) # (i0,c13) -> # of FL associated
    reader = DictReader(open(cluster_csv),delimiter=',')
    for r in reader:
        a,b=r['cluster_id'].split('|')
        cid=b
        pre=a.split('_')[0]
        # current version of SMRTLink does not provide nFL information =____=
        if r['read_type']=='FL':
            fl_ass[(pre,cid)].add(r['read_id'])


    f = open(collapse_prefix + '.abundance.txt', 'w')
    for i in xrange(14): f.write("#\n")
    f.write("pbid\tcount_fl\n")
    for line in open(collapse_prefix + '.group.txt'):
        pbid,members=line.strip().split('\t')
        total = 0
        for m in members.split(','):
            raw = m.split('|')
            cid = raw[1].split('/')[0]
            pre = raw[0].split('_')[0]
            total += len(fl_ass[(pre,cid)])
        f.write("{0}\t{1}\n".format(pbid, total))
    f.close()


def collapse_to_SIRV(out_dir, hq_fastq, cluster_csv, min_count):

    cur_dir = os.getcwd()
    os.chdir(out_dir)

    cmd = "{gmap} -D {gmap_db} -d SIRV -f samse -n 0 -t {cpus} -z sense_force {hq}  > {hq}.sam 2> {hq}.sam.log".format(\
        gmap=GMAP_BIN, gmap_db=GMAP_DB, hq=hq_fastq, cpus=GMAP_CPUS)

    if subprocess.check_call(cmd, shell=True)!=0:
        raise Exception, "ERROR CMD:", cmd

    cmd = "sort -k 3,3 -k 4,4n {hq}.sam > {hq}.sorted.sam".format(hq=hq_fastq)
    if subprocess.check_call(cmd, shell=True)!=0:
        raise Exception, "ERROR CMD:", cmd

    cmd = "collapse_isoforms_by_sam.py --input {hq} --fq -s {hq}.sorted.sam -c 0.99 -i 0.95 \
    --max_fuzzy_junction=5 --dun-merge-5-shorter -o {hq}.no5merge".format(hq=hq_fastq)
    if subprocess.check_call(cmd, shell=True)!=0:
        raise Exception, "ERROR CMD:", cmd

    ### make_abundance_from_CSV
    collapse_prefix = hq_fastq + '.no5merge.collapsed'
    make_abundance_from_Sequel_cluster_csv(cluster_csv, collapse_prefix)



    cmd = "filter_by_count.py {0} --min_count={1}".format(collapse_prefix, min_count)
    if subprocess.check_call(cmd, shell=True)!=0:
        raise Exception, "ERROR CMD:", cmd

    cmd = "filter_away_subset.py {0}.min_fl_{1}".format(collapse_prefix, min_count)
    if subprocess.check_call(cmd, shell=True)!=0:
        raise Exception, "ERROR CMD:", cmd

    os.symlink(collapse_prefix+'.min_fl_'+str(min_count)+'.filtered.abundance.txt', 'touse.count.txt')
    os.symlink(collapse_prefix+'.min_fl_'+str(min_count)+'.filtered.gff', 'touse.gff')
    os.symlink(collapse_prefix+'.group.txt', 'touse.group.txt')

    os.chdir(cur_dir)

def validate_with_SIRV(out_dir, eval_dir):

    out_dir = os.path.abspath(out_dir)
    eval_dir = os.path.abspath(eval_dir)

    cur_dir = os.getcwd()
    os.makedirs(eval_dir)
    os.chdir(eval_dir)
    os.symlink(SIRV_DIR, "SIRV")
    os.symlink(out_dir, "test")

    with open('sample.config', 'w') as f:
        f.write("SAMPLE=SIRV;SIRV/\n")
        f.write("SAMPLE=test;test/\n")
        f.write("\n")
        f.write("GROUP_FILENAME=touse.group.txt\n")
        f.write("GFF_FILENAME=touse.gff\n")
        f.write("COUNT_FILENAME=touse.count.txt\n")

    cmd = "chain_samples.py sample.config count_fl --fuzzy_junction=5"
    if subprocess.check_call(cmd, shell=True)!=0:
        raise Exception, "ERROR CMD:", cmd

    os.chdir(cur_dir)

def eval_result(eval_dir, src_dir, min_count):
    tally = defaultdict(lambda: [])# SIRV --> list of test ids that hit it (can be redundant sometimes due to fuzzy)

    file = os.path.join(eval_dir, 'all_samples.chained_ids.txt')
    assert os.path.exists(file)

    FPs = []
    FNs = []
    for r in DictReader(open(file), delimiter='\t'):
        if r['SIRV'] == 'NA': # is false positive!
            FPs.append(r['test'])
        elif r['test'] == 'NA': # is false negative
            FNs.append(r['SIRV'])
        else:
            tally[r['SIRV']].append(r['test'])

    with open("SIRV_evaluation_summary.txt", 'w') as f:
        f.write("Source Directory: {0}\n".format(src_dir))
        f.write("Evaluation Directory: {0}\n".format(eval_dir))
        f.write("Minimum FL Count cutoff: {0}\n".format(min_count))
        f.write("\n")
        f.write("====================================\n")
        f.write("TP: {0}\n".format(len(tally)))
        f.write("FN: {0}\n".format(len(FNs)))
        f.write("FP: {0}\n".format(len(FPs)))
        f.write("====================================\n")


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser("validate SMRTLink Iso-Seq output against Lexogen SIRV")
    parser.add_argument("smrtlink_dir")
    parser.add_argument("--tmp_dir", default="tmp", help="tmp dirname (default: tmp)")
    parser.add_argument("--eval_dir", default="eval", help="eval dirname (default: eval)")
    parser.add_argument("--min_count", type=int, default=2, help="min FL count cutoff (default:2)")

    args = parser.parse_args()
    o, a, b = link_files(args.smrtlink_dir, args.tmp_dir)
    collapse_to_SIRV(o, a, b, args.min_count)
    validate_with_SIRV(args.tmp_dir, args.eval_dir)
    eval_result(args.eval_dir, args.smrtlink_dir, args.min_count)
