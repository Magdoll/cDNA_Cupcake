import os, sys, gzip, subprocess
import csv
from csv import DictReader, DictWriter
from collections import defaultdict
from multiprocessing import Process

from Bio import SeqIO

import pysam
from bx.intervals.cluster import ClusterTree
import cupcake.io.BioReaders as BioReaders

csv.field_size_limit(100000000)

def collect_cluster_results_multithreaded(group_csv, out_dir, output_prefix, use_BC, chunks):
    indices = set()
    for r in DictReader(open(group_csv), delimiter=','):
        indices.add(r['index'])

    indices = list(indices)
    n = len(indices)
    chunk_size = (n//chunks) + (n%chunks>0)

    pools = []
    for i in range(chunks):
        print("collection worker {0}, starting from {1}".format(i, indices[i*chunk_size]))
        p = Process(target=collect_cluster_results, args=(group_csv,
                                                          out_dir,
                                                          output_prefix+'.'+str(i),
                                                          use_BC,
                                                          indices[(i*chunk_size):(i+1)*chunk_size],))
        p.start()
        pools.append(p)

    for p in pools:
        p.join()

    # fasta file concatenate
    fasta_files = [output_prefix+'.'+str(i)+'.clustered_transcript.rep.fasta' for i in range(chunks)]
    csv_files = [output_prefix+'.'+str(i)+'.clustered_transcript.csv' for i in range(chunks)]

    cmd_cat_fasta = "cat {0} > {1}.clustered_transcript.rep.fasta".format(" ".join(fasta_files), output_prefix)
    cmd_cat_csv1 = "echo \"index,UMI,BC,locus,cluster,ccs_id\" > {0}.clustered_transcript.header".format(output_prefix)
    cmd_cat_csv2 = "cat {0}.clustered_transcript.header {1} > {0}.clustered_transcript.csv".format(output_prefix, " ".join(csv_files))

    try:
        subprocess.check_call(cmd_cat_fasta, shell=True)
    except subprocess.CalledProcessError:
        print("Trouble running CMD: {0}".format(cmd_cat_fasta), file=sys.stderr)
    try:
        subprocess.check_call(cmd_cat_csv1, shell=True)
    except subprocess.CalledProcessError:
        print("Trouble running CMD: {0}".format(cmd_cat_csv1), file=sys.stderr)
    try:
        subprocess.check_call(cmd_cat_csv2, shell=True)
    except subprocess.CalledProcessError:
        print("Trouble running CMD: {0}".format(cmd_cat_csv2), file=sys.stderr)

def collect_cluster_results(group_csv, out_dir, output_prefix, use_BC=False, indices_to_use=None):
    """
    <index>,<UMI>,<locus>,<transcript or NA>,<ccs_id>
    # HQ/LQ transcript IDs go from
    transcript/0 --> 1-TCAAGGTC-transcript/0
    # however singletons (shouldn't have much post cluster) will keep the CCS ID which is ok
    """
    f_out_rep = open(output_prefix + '.clustered_transcript.rep.fasta', 'w')
    f_out_not = open(output_prefix + '.clustered_transcript.not_rep.fasta', 'w')
    f_csv = open(output_prefix + '.clustered_transcript.csv', 'w')
    writer = DictWriter(f_csv, fieldnames=['index', 'UMI', 'BC', 'locus', 'cluster', 'ccs_id'], delimiter=',')
    #writer.writeheader()

    bad = []
    for r in DictReader(open(group_csv), delimiter=','):
        index = r['index']
        if indices_to_use is not None and index not in indices_to_use:
            continue
        umi_key = r['UMI']
        if use_BC:
            umi_key += '-' + r['BC']
        members = list(set(r['members'].split(',')))

        if len(members) == 1: # singleton, no directory, continue
            continue
        d = os.path.join(out_dir, index, umi_key)
        hq = os.path.join(d, 'output.hq.fasta.gz')
        lq = os.path.join(d, 'output.lq.fasta.gz')
        single = os.path.join(d, 'output.singletons.fasta.gz')
        report = os.path.join(d, 'output.cluster_report.csv')

        if not os.path.exists(d):
            print("{0} does not exist! DEBUG mode, ok for now".format(d))
            bad.append(d)
            continue

        if not os.path.exists(report):
            print("{0} does not exist! DEBUG mode, ok for now".format(report))
            bad.append(report)
            continue

        hq_count = 0
        if os.path.exists(hq):
            with gzip.open(hq, 'rt') as handle:
                for seqrec in SeqIO.parse(handle, 'fasta'):
                    if hq_count == 0:
                        f_out_rep.write(">{i}-{u}-{n}\n{s}\n".format(i=index, u=umi_key, n=seqrec.id, s=seqrec.seq))
                    else:
                        f_out_not.write(">{i}-{u}-{n}\n{s}\n".format(i=index, u=umi_key, n=seqrec.id, s=seqrec.seq))
                    hq_count += 1

        if os.path.exists(lq):
            with gzip.open(lq, 'rt') as handle:
                for seqrec in SeqIO.parse(handle, 'fasta'):
                    f_out_not.write(">{i}-{u}-{n}\n{s}\n".format(i=index, u=umi_key, n=seqrec.id, s=seqrec.seq))

        info = {'index': index,
                'UMI': r['UMI'],
                'BC': r['BC'] if use_BC else 'NA',
                'locus': r['locus']}
        for ccs_rec in DictReader(open(report), delimiter=','):
            # cluster_id,read_id,read_type
            # transcript/0,m64012_191109_035807/102828567/ccs,FL
            # transcript/0,m64012_191109_035807/121700628/ccs,FL
            info['cluster'] = ccs_rec['cluster_id']
            info['ccs_id'] = ccs_rec['read_id']
            writer.writerow(info)

        # NOTE: singletons are not part of output.cluster_report.csv
        info['cluster'] = 'NA'
        if os.path.exists(single):
            with gzip.open(single, 'rt') as handle:
                for seqrec in SeqIO.parse(handle, 'fasta'):
                    f_out_not.write(">{i}-{u}-{n}\n{s}\n".format(i=index, u=umi_key, n=seqrec.id, s=seqrec.seq))
                    info['ccs_id'] = seqrec.id
                    writer.writerow(info)

    f_out_not.close()
    f_out_rep.close()
    f_csv.close()
    return bad, f_out_rep.name, f_csv.name

def write_records_to_bam_multithreaded(flnc_tagged_bam, map_seqid_to_group, output_prefix, cpus=1, chunks=4):

    print("Reading {0}...".format(flnc_tagged_bam), file=sys.stdout)
    reader = pysam.AlignmentFile(flnc_tagged_bam, 'rb', check_sq=False)
    read_dict = dict((r.qname, r) for r in reader)

    group_map = defaultdict(lambda: set())   # out_dir_group --> set of seqids
    for seqid,group_name in map_seqid_to_group.items():
        group_map[group_name].add(seqid)

    # now make set back to list
    for k in group_map:
        group_map[k] = list(group_map[k])

    # split group_map into N chunks
    group_map_keys = list(group_map.keys())
    n = len(group_map_keys)
    chunk_size = n//chunks + (n%chunks>0)
    pools = []
    for i in range(chunks):
        p = Process(target=write_records_to_bam, args=(read_dict,
                                                       reader.header,
                                                       group_map,
                                                       output_prefix+'.cmd_'+str(i),
                                                       output_prefix+'.singleton_'+str(i)+'.fasta',
                                                       cpus,
                                                       group_map_keys[(i*chunk_size):((i+1)*chunk_size)],))
        p.start()
        pools.append(p)

    for p in pools:
        p.join()



def write_records_to_bam(read_dict, reader_header, group_map, cmd_filename, singleton_out_fasta, cpus, group_map_keys=None):
    f_cmd = open(cmd_filename, 'w')
    f_singleton = open(singleton_out_fasta, 'w')
    # make the directories + write out flnc.bam

    if group_map_keys is None:
        group_map_keys = group_map.keys()

    single = 0
    for group_name in group_map_keys:
        assert len(group_map[group_name]) == len(set(group_map[group_name])) # DEBUG only can remove later
        single += len(group_map[group_name])==1
    print("{0}-single:{1},not single:{2}".format(cmd_filename,single, len(group_map)-single))
    #return

    for group_name in group_map_keys:
        seqids = group_map[group_name]
        if len(seqids) == 1:  # singleton, just write to the output fasta to save time
            _outdir, _index, _umi_bc = group_name.split('/')
            r = read_dict[seqids[0]]
            # ex: newid 1-TGATACCAC-TGCTTAAAAAAA-m64019_190830_195852/130941517/ccs
            newid = "{i}-{u}-{d}".format(i=_index, u=_umi_bc, d=r.qname)
            f_singleton.write(">{0}\n{1}\n".format(newid, r.seq))
        else:
            print("Writing to {0}/flnc_tagged.bam...".format(group_name), file=sys.stdout)
            os.makedirs(group_name)
            f_out = pysam.AlignmentFile(os.path.join(group_name,"flnc_tagged.bam"), 'wb', header=reader_header)
            for seqid in seqids:
                f_out.write(read_dict[seqid])
            f_out.close()
            f_cmd.write("isoseq3 cluster {d}/flnc_tagged.bam {d}/output.bam --use-qvs --singletons -j {cpus}\n".format(\
                d=group_name, cpus=cpus))
    f_cmd.close()
    f_singleton.close()

def sep_by_UMI(records, umi_bc_dict, out_dir, f_out, map_seqid_to_group, loci_index, use_BC=False):
    """
    :param records: list of SAM records that are overlapping in loci (could be diff strand)
    :param umi_bc_dict: dict of <ccs_id> --> UMI & BC info

    continuously write out to f_out
    <unique_index>, <chr:start-end>, <UMI>, <comma list of ccs_id>

    later we should create <out_dir>/<unique_index>/<UMI>/flnc_tagged.bam  <-- can run isoseq3 cluster on this later
    """
    chrom = records[0].sID
    # first we separate by loci (don't worry about strand)
    tree = ClusterTree(0,0)
    for i,r in enumerate(records):
        tree.insert(r.sStart, r.sEnd, i)

    # now for each region we group by UMI
    for s, e, indices in tree.getregions():
        umi_groups = defaultdict(lambda: set())  # umi(+BC optional) --> set of record indices
        for i in indices:
            if records[i].qID in map_seqid_to_group:  # SPECIAL CASE - split mapping despite using --secondary=no, so ignore the second time we see it
                continue
            if use_BC:
                umi_key = umi_bc_dict[records[i].qID]['UMI']+'-'+umi_bc_dict[records[i].qID]['BC']
            else:
                umi_key = umi_bc_dict[records[i].qID]['UMI']
            umi_groups[umi_key].add(i)  # by using set instead of list, we eliminate chimeric mappings on the same locus

        for umi_key,members_indices in umi_groups.items():
            if use_BC:
                umi, bc = umi_key.split('-')
            else:
                umi, bc = umi_key, 'NA'
            member_ids = list(set([records[i].qID for i in set(members_indices)]))
            info = {'index': loci_index,
                    'UMI': umi,
                    'BC': bc,
                    'locus': "{0}:{1}-{2}".format(chrom,s,e),
                    'size': len(members_indices),
                    'members': ",".join(member_ids)}
            f_out.writerow(info)
            d = os.path.join(out_dir, str(loci_index), umi_key)
            map_seqid_to_group.update(dict((m,d) for m in member_ids))
        loci_index += 1
    return loci_index

def iter_sorted_gmap_record(sam_filename, umi_bc_dict, out_dir, f_out, use_BC=False):
    """
    :param sam_filename: sorted SAM file of tagged FLNC mapped to genome
    :param umi_bc_dict: dict of ccs_id --> dict of UMI/BC/info
    :param out_dir: output directory
    :param f_out: DictWriter object for writing out ccs_id --> group assignment
    :param use_BC: is single cell so also use the "BC" field in addition to "UMI" field
    :return: map_seqid_to_group which is dict of seqid --> group name

    A group is FLNCs that have the same (mapped locus, UMI)
    group name is currently a string of the directory we will create later, which is
       <out_dir>/<loci_index>/<UMI>-<BC>/flnc_tagged.bam
    """
    map_seqid_to_group = {}  # seqid (FLNC CCS id) --> group name
    reader = BioReaders.GMAPSAMReader(sam_filename, True)

    # find first acceptably mapped read
    for r in reader:
        if r.sID != '*': break
    records = [r]
    max_end = r.sEnd

    loci_index = 1
    for r in reader:
        if r.sID == '*': continue
        if r.sID == records[0].sID and r.sStart < records[-1].sStart:
            print("SAM file is NOT sorted. ABORT!", file=sys.stderr)
            sys.exit(-1)
        if r.sID != records[0].sID or r.sStart > max_end:
            print("processing {0}:{1}...{2} records".format(r.sID, max_end, len(records)), file=sys.stdout)
            loci_index = sep_by_UMI(records, umi_bc_dict, out_dir, f_out, map_seqid_to_group, loci_index, use_BC)
            records = [r]
            max_end = r.sEnd
        else:
            records.append(r)
            max_end = max(max_end, r.sEnd)

    return map_seqid_to_group


def run(args):
    print("Reading UMI-BC CSV {0}...".format(args.umi_bc_csv))
    umi_bc_dict = dict((r['id'],r) for r in DictReader(open(args.umi_bc_csv),delimiter='\t'))

    f = open(args.output_prefix + '_founder_info.csv', 'w')
    f_out = DictWriter(f, fieldnames=['index', 'UMI', 'BC', 'locus' ,'size', 'members'])
    f_out.writeheader()
    m = iter_sorted_gmap_record(sam_filename=args.sorted_sam,
                                umi_bc_dict=umi_bc_dict,
                                out_dir=args.out_dir,
                                f_out=f_out,
                                use_BC=args.useBC)

    write_records_to_bam_multithreaded(args.flnc_bam, m, args.output_prefix, cpus=1, chunks=20)
    f.close()

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser("Cluster reads by UMI/BC")
    parser.add_argument("flnc_bam", help="FLNC BAM filename")
    parser.add_argument("sorted_sam", help="Mapped, sorted FLNC SAM filename")
    parser.add_argument("umi_bc_csv", help="Clipped UMI/BC CSV filename")
    parser.add_argument("output_prefix", help="Output prefix")
    parser.add_argument("-d", "--out_dir", default='tmp', help="Cluster out directory (default: tmp/)")
    parser.add_argument("--useBC", default=False, action="store_true", help="Has single cell BC (default: off)")

    args = parser.parse_args()
    run(args)
