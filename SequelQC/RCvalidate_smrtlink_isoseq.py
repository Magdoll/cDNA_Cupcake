__author__ = 'etseng@pacb.com'
__version__ = '1.2'

import os, sys, subprocess
from csv import DictReader
from collections import defaultdict
import SIRVvalidate_smrtlink_isoseq as smrtlink

GENCODE_GTF = '/home/UNIXHOME/etseng/share/gencode/gencode.v29.annotation.gtf'
CAGE_BED = '/home/UNIXHOME/etseng/share/FANTOM/hg38.cage_peak_phase1and2combined_coord.bed'
INTRONPOLIS = '/pbi/dept/bifx/customer_data/Public_Intronpolis/*min_count_10.modified'

SQANTI_QC = '/home/UNIXHOME/etseng/GitHub/SQANTI2/sqanti_qc2.py'
SQANTI_FILTER = '/home/UNIXHOME/etseng/GitHub/SQANTI2/sqanti_filter2.py'

# def make_abundance_from_Sequel_cluster_csv(cluster_csv, collapse_prefix):
#     """
#     cluster_id,read_id,read_type
#     i0_ICE_samplee5686f|c1,m54011_160718_221804/53150541/967_60_CCS,FL
#     """
#     fl_ass = defaultdict(lambda: set()) # (i0,c13) -> # of FL associated
#     all_fl_read_ids = set()
#     reader = DictReader(open(cluster_csv),delimiter=',')
#     for r in reader:
#         a,b=r['cluster_id'].split('|')
#         cid=b
#         pre=a.split('_')[0]
#         # current version of SMRTLink does not provide nFL information =____=
#         if r['read_type']=='FL':
#             fl_ass[(pre,cid)].add(r['read_id'])
#             all_fl_read_ids.add(r['read_id'])
#
#
#     f2 = open(collapse_prefix + '.read_stat.txt', 'w')
#     f2.write("id\tlength\tis_fl\tstat\tpbid\n")
#     f = open(collapse_prefix + '.abundance.txt', 'w')
#     for i in xrange(14): f.write("#\n")
#     f.write("pbid\tcount_fl\n")
#     for line in open(collapse_prefix + '.group.txt'):
#         pbid,members=line.strip().split('\t')
#         total = 0
#         for m in members.split(','):
#             raw = m.split('|')
#             cid = raw[1].split('/')[0]
#             pre = raw[0].split('_')[0]
#             total += len(fl_ass[(pre,cid)])
#             # write to read_stat.txt
#             for fl_read_id in fl_ass[(pre,cid)]:
#                 raw = fl_read_id.split('/')[-1].split('_')
#                 _len = abs(int(raw[0])-int(raw[1]))
#                 f2.write("{0}\t{1}\tY\tunique\t{2}\n".format(fl_read_id, _len, pbid))
#                 all_fl_read_ids.remove(fl_read_id)
#         f.write("{0}\t{1}\n".format(pbid, total))
#     f.close()
#
#     # write the unassigned FL to read_stat.txt
#     for fl_read_id in all_fl_read_ids:
#         raw = fl_read_id.split('/')[-1].split('_')
#         _len = abs(int(raw[0])-int(raw[1]))
#         f2.write("{0}\t{1}\tY\tNA\tNA\n".format(fl_read_id, _len))
#     f2.close()

def collapse_to_hg38(out_dir, hq_fastq, cluster_csv, min_count, aligner_choice, isoseq_version, use_hg19_instead=False):

    cur_dir = os.getcwd()
    os.chdir(out_dir)

    if aligner_choice == 'gmap':
        ref = 'hg19' if use_hg19_instead else 'hg38'
        cmd = "{gmap} -D {gmap_db} -d {ref} -f samse -n 0 -t {cpus} -z sense_force {hq}  > {hq}.sam 2> {hq}.sam.log".format(\
            gmap=smrtlink.GMAP_BIN, gmap_db=smrtlink.GMAP_DB, hq=hq_fastq, cpus=smrtlink.GMAP_CPUS, ref=ref)
    elif aligner_choice == 'minimap2':
        ref = smrtlink.HG19_GENOME if use_hg19_instead else smrtlink.HG38_GENOME
        cmd = "minimap2 -t {cpus} -ax splice -uf --secondary=no -C5  {ref} {hq} > {hq}.sam 2> {hq}.sam.log".format(\
            cpus=smrtlink.GMAP_CPUS, ref=ref, hq=hq_fastq)

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
    cmd = "get_abundance_post_collapse.py " + collapse_prefix + " cluster_report.csv"
    if subprocess.check_call(cmd, shell=True)!=0:
        raise Exception, "ERROR CMD:", cmd

    cmd = "filter_by_count.py {0} --dun_use_group_count --min_count={1}".format(collapse_prefix, min_count)
    if subprocess.check_call(cmd, shell=True)!=0:
        raise Exception, "ERROR CMD:", cmd

    cmd = "filter_away_subset.py {0}.min_fl_{1}".format(collapse_prefix, min_count)
    if subprocess.check_call(cmd, shell=True)!=0:
        raise Exception, "ERROR CMD:", cmd

    rep = collapse_prefix + ".min_fl_{0}.filtered".format(min_count)

    if aligner_choice=='gmap':
        ref = 'hg19' if use_hg19_instead else 'hg38'
        cmd = "{gmap} -D {gmap_db} -d {ref} -f samse -n 0 -t {cpus} -z sense_force {rep}.rep.fq  > {rep}.rep.fq.sam 2> {rep}.rep.fq.sam.log".format(\
            gmap=smrtlink.GMAP_BIN, gmap_db=smrtlink.GMAP_DB, cpus=smrtlink.GMAP_CPUS, rep=rep, ref=ref)
    elif aligner_choice=='minimap2':
        ref = smrtlink.HG19_GENOME if use_hg19_instead else smrtlink.HG38_GENOME
        cmd = "minimap2 -t {cpus} -ax splice -uf --secondary=no -C5 -O6,24 -B4 {ref} {rep}.rep.fq > {rep}.rep.fq.sam 2> {rep}.rep.fq.sam.log".format(\
            cpus=smrtlink.GMAP_CPUS, ref=ref, rep=rep)

    if subprocess.check_call(cmd, shell=True)!=0:
        raise Exception, "ERROR CMD:", cmd

    cmd = "sort -k 3,3 -k 4,4n {rep}.rep.fq.sam > {rep}.rep.fq.sorted.sam".format(rep=rep)
    if subprocess.check_call(cmd, shell=True)!=0:
        raise Exception, "ERROR CMD:", cmd

    os.symlink(rep+'.abundance.txt', 'touse.count.txt')
    os.symlink(rep+'.gff', 'touse.gff')
    os.symlink(collapse_prefix+'.group.txt', 'touse.group.txt')
    os.symlink(rep+'.rep.fq', 'touse.rep.fq')
    os.symlink(rep+'.rep.fq.sam', 'touse.rep.fq.sam')
    os.symlink(rep+'.rep.fq.sorted.sam', 'touse.rep.fq.sorted.sam')

    os.chdir(cur_dir)

def validate_with_Gencode(out_dir, eval_dir, use_hg19_instead=False):
    """
    Run SQANTI to compare with GENCODE
    """

    out_dir = os.path.abspath(out_dir)
    eval_dir = os.path.abspath(eval_dir)

    cur_dir = os.getcwd()
    os.makedirs(eval_dir)
    os.chdir(eval_dir)
    os.symlink(GENCODE_GTF, os.path.basename(GENCODE_GTF))
    os.symlink(os.path.join(out_dir, 'touse.rep.fq.sorted.sam'), 'touse.rep.fq.sorted.sam')
    os.symlink(os.path.join(out_dir, 'touse.rep.fq'), 'touse.rep.fq')

    cmd = "python {sqanti} -t {cpus} {extra} touse.rep.fq {gtf} {genome}".format(\
        sqanti=SQANTI_QC, cpus=smrtlink.GMAP_CPUS, gtf=GENCODE_GTF, \
        genome=smrtlink.HG19_GENOME if use_hg19_instead else smrtlink.HG38_GENOME,
		extra="--cage_peak " + CAGE_BED + " -c " + INTRONPOLIS if not use_hg19_instead else '')

    if subprocess.check_call(cmd, shell=True)!=0:
        raise Exception, "ERROR CMD:", cmd

    cmd = "python {sqanti_filter} touse.rep_classification.txt touse.rep.renamed.fasta touse.rep.renamed_corrected.sam".format(\
        sqanti_filter=SQANTI_FILTER)

    if subprocess.check_call(cmd, shell=True)!=0:
        raise Exception, "ERROR CMD:", cmd

    os.chdir(cur_dir)


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser("validate SMRTLink Iso-Seq output against Lexogen SIRV")
    parser.add_argument("smrtlink_dir")
    parser.add_argument("--tmp_dir", default="tmp", help="tmp dirname (default: tmp)")
    parser.add_argument("--eval_dir", default="eval", help="eval dirname (default: eval)")
    parser.add_argument("--min_count", type=int, default=2, help="min FL count cutoff (default:2)")
    parser.add_argument("--aligner_choice", default='minimap2', choices=('gmap', 'star', 'minimap2'), help="Aligner choice (default: minimap2)")
    parser.add_argument("--use_hg19_instead", default=False, action="store_true")

    args = parser.parse_args()
    o, a, b, v = smrtlink.link_files(args.smrtlink_dir, args.tmp_dir)
    collapse_to_hg38(o, a, b, args.min_count, args.aligner_choice, v, args.use_hg19_instead)
    validate_with_Gencode(args.tmp_dir, args.eval_dir, args.use_hg19_instead)
