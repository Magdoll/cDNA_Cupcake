__author__ = 'etseng@pacb.com'
__version__ = '2.0'

import os, sys, subprocess, shutil
from csv import DictReader
from collections import defaultdict

import SIRVvalidate_smrtlink_isoseq as smrtlink_sirv

"""
Files that should already be in the same directory:

--- hq_isoforms.fasta
--- lq_isoforms.fasta
--- cluster_report.csv


"""

GMAP_BIN = smrtlink_sirv.GMAP_BIN
GMAP_DB = smrtlink_sirv.GMAP_DB
GMAP_CPUS = smrtlink_sirv.GMAP_CPUS
SIRV_DIR = smrtlink_sirv.SIRV_DIR
SIRV_GENOME = smrtlink_sirv.SIRV_GENOME

HG38_GENOME = smrtlink_sirv.HG38_GENOME

STAR_BIN = smrtlink_sirv.STAR_BIN
SIRV_STAR_DB = smrtlink_sirv.SIRV_STAR_DB


def link_files(src_dir, out_dir):
    hq_fastq = os.path.join(src_dir, 'hq_isoforms.fastq')
    cluster_csv = os.path.join(src_dir, 'cluster_report.csv')

    o_fq = os.path.join(out_dir, 'hq_isoforms.fastq')
    o_csv = os.path.join(out_dir, 'cluster_report.csv')
    if os.path.exists(out_dir):
        assert len(os.popen("diff {0} {1}".format(hq_fastq, o_fq)).read().strip()) == 0
        print >> sys.stderr, "Re-using data in {0}....".format(out_dir)
    else:
        os.makedirs(out_dir)
        os.symlink(hq_fastq, o_fq)
        os.symlink(cluster_csv, o_csv)

    return out_dir, 'hq_isoforms.fastq', 'cluster_report.csv'


def make_abundance_from_Sequel_cluster_csv(cluster_csv, collapse_prefix):
    """
    cluster_id,read_id,read_type
    cb10265_c1,m54086_170204_081430/7143755/4015_60_CCS,FL
    cb10265_c1,m54086_170204_081430/50069845/2446_53_CCS,NonFL
    cb10265_c1,m54086_170204_081430/69862364/30_8329_CCS,NonFL
    """
    fl_ass = defaultdict(lambda: set()) # cb10265_c1 -> # of FL associated
    nfl_ass = defaultdict(lambda: set()) # cb10265_c1 -> # of nFL associated
    nfl_hit_count = defaultdict(lambda: 0) # nFL seqid --> # of cids it belongs to
    reader = DictReader(open(cluster_csv),delimiter=',')
    for r in reader:
        cid = r['cluster_id']
        # current version of SMRTLink does not provide nFL information =____=
        if r['read_type']=='FL':
            fl_ass[cid].add(r['read_id'])
        else:
            nfl_ass[cid].add(r['read_id'])
            nfl_hit_count[r['read_id']] += 1

    f = open(collapse_prefix + '.abundance.txt', 'w')
    for i in xrange(14): f.write("#\n")
    f.write("pbid\tcount_fl\tcount_nfl\n")
    for line in open(collapse_prefix + '.group.txt'):
        pbid,members=line.strip().split('\t')
        total_fl = 0
        total_nfl = 0
        for m in members.split(','):
            i = m.find('|')
            if i > 0: m = m[i+1:]
            cid = m.split('/')[0]
            total_fl += len(fl_ass[cid])
            total_nfl += len(fl_ass[cid]) + len(filter(lambda x: nfl_hit_count[x]==1, nfl_ass[cid]))
        f.write("{0}\t{1}\t{2}\n".format(pbid, total_fl, total_nfl))
    f.close()

def sanity_check_script_dependencies():
    if os.system("chain_samples.py -h > /dev/null")!=0:
        print >> sys.stderr, "chain_samples.py required in PATH! Please install Cupcake ToFU!"
        sys.exit(-1)

def collapse_to_SIRV(out_dir, hq_fastq, cluster_csv, min_count, aligner_choice):

    cur_dir = os.getcwd()
    os.chdir(out_dir)

    # don't re-do alignment if already there
    if os.path.exists(hq_fastq+'.sam'):
        os.remove('touse.group.txt')
        os.remove('touse.count.txt')
        os.remove('touse.gff')
    else:
        if aligner_choice=='gmap':
            cmd = "{gmap} -D {gmap_db} -d SIRV -f samse -n 0 -t {cpus} -z sense_force {hq}  > {hq}.sam 2> {hq}.sam.log".format(\
            gmap=GMAP_BIN, gmap_db=GMAP_DB, hq=hq_fastq, cpus=GMAP_CPUS)
        elif aligner_choice=='star':
            cmd = "{star} {db} {hq} {hq}.sam --cpus {cpus}".format(\
                star=STAR_BIN, db=SIRV_STAR_DB, hq=hq_fastq, cpus=GMAP_CPUS)
        elif aligner_choice=='minimap2':
            cmd = "minimap2 -ax splice -uf --secondary=no --splice-flank=no -C5 -t {cpus} {ref} {hq} > {hq}.sam 2> {hq}.sam.log".format(\
                cpus=GMAP_CPUS, ref=SIRV_GENOME, hq=hq_fastq)
        else:
            raise Exception, "Unrecognized aligner choice: {0}!".format(aligner_choice)

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
    if os.path.exists(eval_dir):
        shutil.rmtree(eval_dir)
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
    parser.add_argument("--aligner_choice", default='star', choices=('gmap', 'star', 'minimap2'), help="Aligner choice (default: star)")

    args = parser.parse_args()
    sanity_check_script_dependencies()

    o, a, b = link_files(os.path.realpath('.'), args.tmp_dir)
    collapse_to_SIRV(o, a, b, args.min_count, args.aligner_choice)
    validate_with_SIRV(args.tmp_dir, args.eval_dir)
    eval_result(args.eval_dir, args.smrtlink_dir, args.min_count)
