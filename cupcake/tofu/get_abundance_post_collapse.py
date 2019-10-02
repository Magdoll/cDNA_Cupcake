__author__ = 'etseng@pacb.com'

"""
After running collapse, combine with the cluster FL/nFL assignment to get abundance information.
Compliant with SMRTLink/SMRTAnalysis 3.2+ formatting.
May break on older SMRTPortal/SMRTAnalysis 2.x formatting.

cluster_report.csv format in IsoSeq1:

cluster_id,read_id,read_type (FL or NonFL)
i0_ICE_sample03146a|c4,m54026_161015_230744/74318463/31_970_CCS,FL

NOTE: in cluster_report, the sample prefix is i<bin>_ICE_<sample>, but in HQ isoforms is i<bin>_HQ_<sample>!
Must address this :P

cluster_report.csv format in IsoSeq2:

cluster_id,read_id,read_type
cb2729_c2,m54056_171130_193019/8388911/26_1296_CCS,FL


cluster_report.csv format in IsoSeq3:

cluster_id,read_id,read_type
transcript/0,m54056_171130_193019/8388911/ccs,FL


where .group.txt is:

PBfusion.1      HQ_sample0ZPg9hS7|cb7607_c93041/f2p0/1713
PBfusion.2      HQ_sample0ZPg9hS7|cb7607_c16635/f3p0/810,HQ_sample0ZPg9hS7|cb7607_c32934/f2p1/1066

"""


import os, sys, re
from collections import defaultdict
from csv import DictReader

def get_roi_len(seqid):
    # before isoseq3: <movie>/<zmw>/<start>_<end>_CCS
    # for isoseq3: <movie>/<zmw>/ccs
    if seqid.endswith('/ccs'):
        print("WARNING: isoseq3 format detected. Output `length` column will be `NA`.", file=sys.stderr)
        return "NA"
    elif not seqid.endswith('_CCS'):
        print("Sequence ID format must be <movie>/<zmw>/<start>_<end>_CCS or <movie>/<zmw>/ccs! Abort!", file=sys.stderr)
        sys.exit(-1)
    s, e, junk = seqid.split('/')[2].split('_')
    return abs(int(s)-int(e))



cluster_rex_sa3 = re.compile('(i\d+_ICE_\S+\|c\d+)')
cluster_rex_sa2 = re.compile('(i\d+[HL]Q_\S+|c\d+)')
# ex: cb100_c7,m54006_170206_215027/70189550/2661_56_CCS,FL
cluster_rex_tofu2 = re.compile('(cb\d+_c\d+)')
cluster_rex_isoseq3 = re.compile('(transcript/\d+)') # for isoseq3: transcript/0, transcript/1, etc
cluster_rex_isoseq3_mapped = re.compile('\S+(transcript/\d+)') # for isoseq post SL mapping


def read_group_filename(group_filename, is_cid=True):
    """
    Make the connection between partitioned results and final (ex: PB.1.1)
    The partitioned results could either be ICE cluster (ex: i1_c123) or RoIs

    Group filename format: (SMRTLink/SA3.x, IsoSeq1)
    PB.1.1  i0_HQ_sample03146a|c18072/f2p30/279,i0_HQ_sample03146a|c38435/f2p0/180
    PB.2.1  i2_HQ_sample03146a|c46363/f22p0/2776
    PB.2.2  i2_HQ_sample03146a|c46572/f37p0/2761

    For IsoSeq2:
    PBfusion.1      HQ_sample0ZPg9hS7|cb7607_c93041/f2p0/1713
    PBfusion.2      HQ_sample0ZPg9hS7|cb7607_c16635/f3p0/810,HQ_sample0ZPg9hS7|cb7607_c32934/f2p1/1066

	For IsoSeq3, two possible flavors depending on using conda version or SL Mapping version:
	PB.1.1	transcript/0,transcript/1
	(or)
	PB.1.1	Sample1_RC0_50pM_transcript/0,Sample1_RC0_50pM_transcript/1   <-- in this case we want to strip the prefix out, just keep 'transcript/0'

    Note: LQ/HQ isoforms use the _HQ_ or _LQ_ prefix, but in cluster_report.csv it will be _ICE_
    So we just replace _HQ_ or _LQ_ --> both to be _ICE_

    Return: dict of seq_or_ice_cluster --> collapsed cluster ID
          (ex: in IsoSeq3 it is 'transcript/0' --> 'PB.1.1')
    """
    cid_info = {} # ex: i1 --> c123 --> PB.1.1, or c123 --> PB.1.1

    for line in open(group_filename):
        pbid, members = line.strip().split('\t')
        for cid in members.split(','):
            m = cluster_rex_isoseq3_mapped.match(cid)
            if m is not None: # IsoSeq3 after mapping, strip it down to 'transcript/0'
                cid = m.group(1)
            else:
                m = cluster_rex_isoseq3.match(cid)
                if m is not None:
                    pass # pass, nothing to do since the cid is correctly 'transcript/0'
                else:
                    if is_cid: cid = cid.split('/')[0]
                    m = cluster_rex_tofu2.match(cid)
                    if m is not None:
                        cid = m.group(1)
                    else:
                        m = cluster_rex_sa2.match(cid)
                        if m is not None:
                            cid = m.group(1)
                            # replace _HQ_ or _LQ_ with _ICE_ to be compatible with cluster_report.csv later
                            cid = cid.replace('_HQ_', '_ICE_')
                            cid = cid.replace('_LQ_', '_ICE_')
                        elif cluster_rex_sa3.match(cid) is not None:
                            pass # nothing to do, ID is good
                        elif cid.startswith('HQ_') or cid.startswith('LQ_'): # isoseq2
                            cid = cid.split('|')[1].split('/')[0]  # cid = cb7607_c93041, for example
                        else:
                            raise Exception("Unrecognized id format {0} in {1}!".format(cid, group_filename))
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
        f.write("id\tlength\tis_fl\tstat\tpbid\n")
    elif output_mode == 'a':
        f = open(output_filename, 'a')
    else:
        raise Exception("Output mode {0} not valid!".format(output_mode))

    unmapped_holder = set()
    for r in DictReader(open(csv_filename), delimiter=','):
        cid = str(r['cluster_id'])
        m = cluster_rex_sa3.match(cid)
        if m is None:
            m = cluster_rex_sa2.match(cid)
            if m is None:
                m = cluster_rex_tofu2.match(cid)
                if m is None:
                    m = cluster_rex_isoseq3.match(cid)
                    if m is not None:
                        cid = m.group(1) # make sure cid is transcript/0, transcript/1, without any prefix before 'transcript'
                    else:
                        raise Exception("cluster_id {0} is not a valid cluster ID!".format(cid))

        x = r['read_id']
        if cid in cid_info:
            # if is FL, must be unique
            if r['read_type'] == 'FL':
                pbid, stat = cid_info[cid], 'unique'
                f.write("{id}\t{len}\t{is_fl}\t{stat}\t{pbid}\n".format(\
                    id=x, len=get_roi_len(x), is_fl='Y', stat=stat, pbid=pbid))
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
                f.write("{id}\t{len}\t{is_fl}\t{stat}\t{pbid}\n".format(\
                    id=x, len=get_roi_len(x), is_fl='Y', stat=stat, pbid=pbid))
            else: # nonFL could be multi-mapped, so not sure it is truly unmapped, put in holder
                unmapped_holder.add(x)

    # now we can go through the list of mapped/unmapped to see which are uniquely mapped which are not
    for seqid, pbids in mapped.items():
        if len(pbids) == 1: # unique
            stat = 'unique'
        else:
            stat = 'ambiguous'
        for pbid in pbids:
            f.write("{id}\t{len}\t{is_fl}\t{stat}\t{pbid}\n".format(
                id=seqid, len=get_roi_len(seqid), is_fl='N', stat=stat, pbid=pbid))

    unmapped_holder = unmapped_holder.difference(mapped)
    pbid, stat = 'NA', 'unmapped'
    for x in unmapped_holder:
        f.write("{id}\t{len}\t{is_fl}\t{stat}\t{pbid}\n".format(
            id=x, len=get_roi_len(x), is_fl='N', stat=stat, pbid=pbid))

    f.close()


def make_abundance_file(read_count_filename, output_filename, given_total=None, restricted_movies=None, write_header_comments=True):
    """
    If given_total is not None, use it instead of the total count based on <read_count_filename>
    given_total should be dict of {fl, nfl, nfl_amb}

    total_ids['fl'] --- number of FL, mapped or unmapped
    total_ids['nfl'] --- number of mapped unique nFL plus unmapped nFL, but excludes mapped ambiguous nFL
    total_ids['nfl_amb'] --- number of mapped, non-unique (ambiguous) nFL
    """
    total_ids = {'fl': set(), 'nfl': set(), 'nfl_amb': set()}
    tally = defaultdict(lambda: {'fl': 0, 'nfl': 0, 'nfl_amb': 0}) # pbid, could be NA --> # of FL reads mapped to it
    amb_count = defaultdict(lambda: []) # non-fl id --> list of pbid matches

    reader = DictReader(open(read_count_filename), delimiter='\t')
    for r in reader:
        movie = r['id'].split('/')[0]
        if restricted_movies is None or movie in restricted_movies:
            if r['pbid'] != 'NA':
                if r['is_fl'] == 'Y': # FL, must be uniquely mapped
                    assert r['stat'] == 'unique'
                    tally[r['pbid']]['fl'] += 1
                    total_ids['fl'].add(r['id'])
                else: # non-FL, can be ambiguously mapped
                    if r['stat'] == 'unique':
                        tally[r['pbid']]['nfl'] += 1
                        total_ids['nfl'].add(r['id'])
                    else:
                        assert r['stat'] == 'ambiguous'
                        amb_count[r['id']].append(r['pbid'])
                        total_ids['nfl_amb'].add(r['id'])
            else: # even if it is unmapped it still counts in the abundance total!
                if r['is_fl'] == 'Y':
                    total_ids['fl'].add(r['id'])
                else:
                    total_ids['nfl'].add(r['id'])


    # put the ambiguous back in tally weighted
    for seqid, pbids in amb_count.items():
        weight = 1. / len(pbids)
        for pbid in pbids:
             tally[pbid]['nfl_amb'] += weight

    if given_total is not None:
        use_total_fl = given_total['fl']
        use_total_nfl = given_total['fl'] + given_total['nfl']
        # ToDo: the below is NOT EXACTLY CORRECT!! Fix later!
        use_total_nfl_amb = given_total['fl'] + given_total['nfl'] + given_total['nfl_amb']
    else:
        use_total_fl = len(total_ids['fl'])
        use_total_nfl = len(total_ids['fl']) + len(total_ids['nfl'])
        use_total_nfl_amb = len(total_ids['fl']) + len(total_ids['nfl']) + len(total_ids['nfl_amb'])

    f = open(output_filename,'w')
    if write_header_comments:
        f.write("#\n")
        f.write("# -----------------\n")
        f.write("# Field explanation\n")
        f.write("# -----------------\n")
        f.write("# count_fl: Number of associated FL reads\n")
        f.write("# count_nfl: Number of associated FL + unique nFL reads\n")
        f.write("# count_nfl_amb: Number of associated FL + unique nFL + weighted ambiguous nFL reads\n")
        f.write("# norm_fl: count_fl / total number of FL reads, mapped or unmapped\n")
        f.write("# norm_nfl: count_nfl / total number of mapped unique or unmapped (FL + nFL) reads\n")
        f.write("# norm_nfl_amb: count_nfl_amb / total number of all reads\n")
        f.write("# Total Number of FL reads: {0}\n".format(use_total_fl))
        f.write("# Total Number of FL + unique nFL reads: {0}\n".format(use_total_nfl))
        f.write("# Total Number of all reads: {0}\n".format(use_total_nfl_amb))
        f.write("#\n")
    f.write("pbid\tcount_fl\tcount_nfl\tcount_nfl_amb\tnorm_fl\tnorm_nfl\tnorm_nfl_amb\n")

    keys = list(tally.keys())
    keys.sort(key=lambda x: list(map(int, x.split('.')[1:]))) # sort by PB.1, PB.2....
    for k in keys:
        v = tally[k]
        a = v['fl']
        b = a + v['nfl']
        c = b + v['nfl_amb']
        f.write("{0}\t{1}\t{2}\t{3}\t".format(k, a, b, c))
        f.write("{0:.4e}\t{1:.4e}\t{2:.4e}\n".format(a*1./use_total_fl, b*1./use_total_nfl, c*1./use_total_nfl_amb))
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

    cid_info = read_group_filename(collapse_prefix + ".group.txt", is_cid=True)

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
