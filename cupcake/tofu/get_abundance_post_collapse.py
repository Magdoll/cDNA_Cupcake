#!/usr/bin/env python3
__author__ = 'etseng@pacb.com'

"""
After running collapse, combine with the cluster FL/nFL assignment to get abundance information.

cluster_report.csv format in IsoSeq3:

cluster_id,read_id,read_type
transcript/0,m54056_171130_193019/8388911/ccs,FL

where .group.txt is:

PB.1.1     transcript/0,transcript1
PB.1.2     transcript/2
"""


import os, sys, re
from collections import Counter
from csv import DictReader, DictWriter


def read_group_filename(group_filename):
    """
    Make the connection between partitioned results and final (ex: PB.1.1)
    Return: dict of seq_or_ice_cluster --> collapsed cluster ID
          (ex: in IsoSeq3 it is 'transcript/0' --> 'PB.1.1')
    """
    cid_info = {} # ex: i1 --> c123 --> PB.1.1, or c123 --> PB.1.1

    for line in open(group_filename):
        pbid, members = line.strip().split('\t')
        for cid in members.split(','):
            cid_info[cid] = pbid

    return cid_info

def output_read_count_IsoSeq_csv(cid_info, csv_filename, output_filename, output_mode='w'):
    """
    Turn into read_stats.txt format:
    id \t length \t is_fl \t stat \t pbid
    """
    mapped = {}  # nFL seq -> list of (sample_prefix, cluster) it belongs to

    if output_mode == 'w':
        f = open(output_filename, 'w')
        f.write("id\tis_fl\tstat\tpbid\n")
    elif output_mode == 'a':
        f = open(output_filename, 'a')
    else:
        raise Exception("Output mode {0} not valid!".format(output_mode))

    unmapped_holder = set()
    for r in DictReader(open(csv_filename), delimiter=','):
        cid = str(r['cluster_id'])
        x = r['read_id']
        if cid in cid_info:
            # if is FL, must be unique
            if r['read_type'] == 'FL':
                pbid, stat = cid_info[cid], 'unique'
                f.write("{id}\t{is_fl}\t{stat}\t{pbid}\n".format(\
                    id=x, is_fl='Y', stat=stat, pbid=pbid))
            else: # nonFL could be multi-mapped, must wait and see
                assert r['read_type'] == 'NonFL'
                # is only potentially unmapped, add all (movie-restricted) members to unmapped holder
                pbid = cid_info[cid]
                if x not in mapped: mapped[x] = set()
                mapped[x].add(pbid)
        else:
            # unmapped, could be either FL or nFL
            # if FL, can immediately write out since it only appears once in cluster_report.csv
            if r['read_type'] == 'FL':
                pbid, stat = 'NA', 'unmapped'
                f.write("{id}\t{is_fl}\t{stat}\t{pbid}\n".format(\
                    id=x, is_fl='Y', stat=stat, pbid=pbid))
            else: # nonFL could be multi-mapped, so not sure it is truly unmapped, put in holder
                unmapped_holder.add(x)

    # now we can go through the list of mapped/unmapped to see which are uniquely mapped which are not
    for seqid, pbids in mapped.items():
        if len(pbids) == 1: # unique
            stat = 'unique'
        else:
            stat = 'ambiguous'
        for pbid in pbids:
            f.write("{id}\t{is_fl}\t{stat}\t{pbid}\n".format(
                id=seqid, is_fl='N', stat=stat, pbid=pbid))

    unmapped_holder = unmapped_holder.difference(mapped)
    pbid, stat = 'NA', 'unmapped'
    for x in unmapped_holder:
        f.write("{id}\\t{is_fl}\t{stat}\t{pbid}\n".format(
            id=x, is_fl='N', stat=stat, pbid=pbid))

    f.close()


def make_abundance_file(read_count_filename, output_filename, given_total=None, restricted_movies=None, write_header_comments=True):
    """
    If given_total is not None, use it instead of the total count based on <read_count_filename>
    """
    tally = Counter() # pbid --> FL count

    reader = DictReader(open(read_count_filename), delimiter='\t')
    for r in reader:
        movie = r['id'].split('/')[0]
        if restricted_movies is None or movie in restricted_movies:
            if r['is_fl']!='Y':
                print("sequence {0} has `is_fl` field set to `{1}`. Ignoring this read.".format(r['id'], r['is_fl']), file=sys.stderr)
                continue
            if r['pbid']=='NA':
                print("sequence {0} has `pbid` field set to `{1}`. Ignoring this read.".format(r['id'], r['pbid']), file=sys.stderr)
                continue
            assert r['stat'] == 'unique'
            tally[r['pbid']] += 1

    if given_total is not None:
        use_total_fl = given_total
    else:
        use_total_fl = sum(tally.values())

    COUNT_FIELDS = ['pbid', 'count_fl', 'norm_fl']
    f = open(output_filename,'w')
    if write_header_comments:
        f.write("#\n")
        f.write("# -----------------\n")
        f.write("# Field explanation\n")
        f.write("# -----------------\n")
        f.write("# count_fl: Number of associated FL reads\n")
        f.write("# norm_fl: count_fl / total number of FL reads, mapped or unmapped\n")
        f.write("# Total Number of FL reads: {0}\n".format(use_total_fl))
        f.write("#\n")

    writer = DictWriter(f, COUNT_FIELDS, delimiter='\t')
    writer.writeheader()

    keys = list(tally.keys())
    keys.sort(key=lambda x: list(map(int, x.split('.')[1:]))) # sort by PB.1, PB.2....
    for k in keys:
        count_fl = tally[k]
        norm_fl  = count_fl*1./use_total_fl
        rec = {'pbid': k, 'count_fl': count_fl, 'norm_fl': "{0:.4e}".format(norm_fl)}
        writer.writerow(rec)
    f.close()




def get_abundance_post_collapse(collapse_prefix, cluster_report_csv, output_prefix, restricted_movies=None):
    """

    :param collapse_prefix: collapse prefix filename (must have .group.txt present)
    :param prefix_dict:
    :param output_prefix:
    :param restricted_movies:
    :return:
    """
    group_filename = collapse_prefix + ".group.txt"
    if not os.path.exists(group_filename):
        print("File {0} does not exist. Abort!".format(group_filename), file=sys.stderr)
        sys.exit(-1)

    if not os.path.exists(cluster_report_csv):
        print("File {0} does not exist. Abort!".format(cluster_report_csv), file=sys.stderr)
        sys.exit(-1)

    cid_info = read_group_filename(collapse_prefix + ".group.txt")

    output_read_count_IsoSeq_csv(cid_info, cluster_report_csv, output_prefix + '.read_stat.txt')
    print("Read stat file written to", output_prefix + '.read_stat.txt', file=sys.stderr)
    make_abundance_file(output_prefix + '.read_stat.txt', output_prefix + '.abundance.txt', restricted_movies=restricted_movies)
    print("Abundance file written to", output_prefix + '.abundance.txt', file=sys.stderr)


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser("Get abundance/read stat information after running collapse script. Works for Iso-Seq1, 2, and 3 output.")
    parser.add_argument("collapse_prefix", help="Collapse prefix (must have .group.txt)")
    parser.add_argument("cluster_report", help="Cluster CSV report")

    args = parser.parse_args()
    get_abundance_post_collapse(args.collapse_prefix, args.cluster_report, args.collapse_prefix)
