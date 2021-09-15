#!/usr/bin/env python
__author__ = 'etseng@pacb.com'

import pdb
import os, sys
import itertools
from pickle import *
from collections import defaultdict, namedtuple

from Bio import SeqIO
from Bio.SeqIO import SeqRecord
from Bio.Seq import Seq

from bx.intervals.cluster import ClusterTree
from cupcake.io import BioReaders
from cupcake.tofu.utils import check_ids_unique
from cupcake.io.SeqReaders import LazyFastaReader, LazyFastqReader
from cupcake.tofu.branch import branch_simple2
from cupcake.tofu.compare_junctions import compare_junctions
from cupcake.tofu.get_abundance_post_collapse import get_abundance_post_collapse


def get_isoform_index(in_order_ranges, chrom, start, end):
    """
    in_order_ranges: list of in order (chrom, start, end)

    NOTE: ranges can shift aftering merging multiple fusion candidates, so we're not asking for precise, just overlap
    """
    for i, (chrom2, start2, end2) in enumerate(in_order_ranges):
        if chrom2==chrom and (start2 <= start < end2 or start <= start2 < end or \
                start2 <= end < end2 or start <= end2 < end):
            return i
    return None


def pick_rep(fa_fq_filename, sam_filename, gff_filename, group_filename, output_filename, fusion_candidate_ranges, is_fq=False):
    """
    For each group, select the representative record
    Always pick the longest one!
    """
    if is_fq:
        fd = LazyFastqReader(fa_fq_filename)
        fout = open(output_filename, 'w')
    else:
        fd = LazyFastaReader(fa_fq_filename)
        fout = open(output_filename, 'w')

    rep_info = {}
    id_to_rep = {}
    for line in open(group_filename):
        pb_id, members = line.strip().split('\t')
        print("Picking representative sequence for", pb_id, file=sys.stdout)
        best_id = None
        best_seq = None
        max_len = 0
        for x in members.split(','):
            if len(fd[x].seq) >= max_len:
                best_id = x
                best_seq = fd[x].seq
                best_qual = fd[x].letter_annotations['phred_quality'] if is_fq else None
                max_len = len(fd[x].seq)
        rep_info[pb_id] = (best_id, best_seq, best_qual)
        id_to_rep[best_id] = pb_id

    f_gff = open(gff_filename, 'w')
    coords = {}
    record_storage = {} # temporary storage for the .1 record to write in conjunction with second record
    for r in BioReaders.GMAPSAMReader(sam_filename, True):
        if r.qID in id_to_rep:
            pb_id = id_to_rep[r.qID]
            # make coordinates & write the SAM file
            isoform_index = get_isoform_index(fusion_candidate_ranges[r.qID], r.sID, r.sStart, r.sEnd)
            if r.qID not in coords:
                coords[r.qID] = [None] * len(fusion_candidate_ranges[r.qID])
                record_storage[pb_id] = [None] * len(fusion_candidate_ranges[r.qID])
            coords[r.qID][isoform_index] = "{0}:{1}-{2}({3})".format(r.sID, r.sStart, r.sEnd, r.flag.strand)
            record_storage[pb_id][isoform_index] = r

    for pb_id, records in record_storage.items():
        for i,r in enumerate(records):
            isoform_index = i + 1
            f_gff.write("{chr}\tPacBio\ttranscript\t{s}\t{e}\t.\t{strand}\t.\tgene_id \"{pi}\"; transcript_id \"{pi}.{j}\";\n".format(\
                chr=r.sID, s=r.segments[0].start+1, e=r.segments[-1].end, pi=pb_id, j=isoform_index, strand=r.flag.strand))
            for s in r.segments:
                f_gff.write("{chr}\tPacBio\texon\t{s}\t{e}\t.\t{strand}\t.\tgene_id \"{pi}\"; transcript_id \"{pi}.{j}\";\n".format(\
                    chr=r.sID, s=s.start+1, e=s.end, pi=pb_id, j=isoform_index, strand=r.flag.strand))
    f_gff.close()

    for pb_id in rep_info:
        best_id, best_seq, best_qual = rep_info[pb_id]
        _id_ = "{0}|{1}|{2}".format(pb_id, "+".join(coords[best_id]), best_id)
        _seq_ = best_seq
        if is_fq:
            SeqIO.write(SeqRecord(_seq_, id=_id_, letter_annotations={'phred_quality':best_qual}), fout, 'fastq')
        else:
            SeqIO.write(SeqRecord(_seq_, id=_id_), fout, 'fasta')

def sep_by_strand(records):
    output = {'+':[], '-':[]}
    for r in records:
        output[r.flag.strand].append(r)
    return output

def is_fusion_compatible(r1, r2, max_fusion_point_dist, max_exon_end_dist, allow_extra_5_exons):
    """
    Helper function for: merge_fusion_exons()

    Check that:
    (1) r1, r2 and both in the 5', or both in the 3'
    (2) if single-exon, fusion point must be close by
        if multi-exon, every junction identical (plus below is True)
    (3) if allow_extra_5_exons is False, num exons must be the same
        if allow_extra_5_exons is True, only allow additional 5' exons
    """
#    _ids = 'i1a_c1603/f67p459/1248,i1b_c19881/f7p368/1235,newClontech_i0HQ|c18279/f6p24/1229,i2b_c22046/f2p494/2157,i2a_c4714/f10p554/2152'.split(',')
#    if r1.qID in _ids or r2.qID in _ids:
#        pdb.set_trace()
    # first need to figure out ends
    # also check that both are in the 5' portion of r1 and r2
    assert r1.flag.strand == r2.flag.strand
    if r1.qStart <= .5*r1.qLen: # in the 5' portion of r1
        if r2.qStart > .5*r2.qLen: # in the 3' portion, reject
            return False
        in_5_portion = True
    else: # in the 3' portion of r1
        if r2.qStart <= .5*r2.qLen:
            return False
        in_5_portion = False
    plus_is_5end = (r1.flag.strand == '+')

    r1.strand = r1.flag.strand
    r2.strand = r2.flag.strand
    type = compare_junctions(r1, r2)
    if type == 'exact':
        if len(r1.segments) == 1:
            if len(r2.segments) == 1:
                # single exon case, check fusion point is close enough
                if in_5_portion and plus_is_5end: dist = abs(r1.sStart - r2.sStart)
                else: dist = abs(r1.sEnd - r2.sEnd)
                return dist <= max_fusion_point_dist
            else:
                raise Exception("Not possible case for multi-exon transcript and " + \
                        "single-exon transcript to be exact!")
        else: # multi-exon case, must be OK
            return True
    if type == 'subset':
            r1, r2 = r2, r1 #  rotate so r1 is always the longer one
    if type == 'super' or type == 'subset':
        n2 = len(r2.segments)
        if allow_extra_5_exons:
            # check that the 3' junction is identical
            # also check that the 3' end is relatively close
            # AND the 5' start of r2 is sandwiched between the matching r1 exon coordinates
            if in_5_portion and plus_is_5end:
                if abs(r1.segments[-1].start - r2.segments[-1].start) > max_exon_end_dist: return False
                if abs(r1.segments[-1].end - r2.segments[-1].end) > max_fusion_point_dist: return False
                return r1.segments[-n2].start <= r2.segments[0].start < r1.segments[-n2].end
            elif in_5_portion and (not plus_is_5end):
                if abs(r1.segments[0].end - r2.segments[0].end) > max_exon_end_dist: return False
                if abs(r1.segments[0].start - r2.segments[0].start) > max_fusion_point_dist: return False
                return r1.segments[n2-1].start <= r2.segments[-1].end < r1.segments[n2-1].end
            else:
                return False
        else: # not OK because number of exons must be the same
            return False
    else: #ex: partial, nomatch, etc...
        return False

def merge_fusion_exons(records, max_fusion_point_dist, max_exon_end_dist, allow_extra_5_exons):
    """
    Records is a list of overlapping GMAP SAM Records (must be on same strand)
    Unlike regular (non-fusion) mapping, only merge records if:

    (1) for multi-exon, every junction is identical
        for single-exon, the fusion point is no bigger than <max_fusion_point_dist> apart

    (2) if allow_extra_5_exons is False, number of exons must be the same
        if allow_extra_5_exons is True, only merge if the extension is in the 5' direction

    Returns a list of grouped records, ex: [[r1,r2], [r3], [r4, r5, r6]]....
    which can be sent to BranchSimple.process_records for writing out
    """
    output = [[records[0]]]
    for r1 in records[1:]:
        merged = False
        # go through output, seeing if mergeable
        for i, r2s in enumerate(output):
#            _ids = 'i1a_c1603/f67p459/1248,i1b_c19881/f7p368/1235,newClontech_i0HQ|c18279/f6p24/1229,i2b_c22046/f2p494/2157,i2a_c4714/f10p554/2152'.split(',')
#            if r1.qID in _ids and any(r2.qID in _ids for r2 in r2s):
#                pdb.set_trace()
            if all(is_fusion_compatible(r1, r2, max_fusion_point_dist, max_exon_end_dist, allow_extra_5_exons) for r2 in r2s):
                output[i].append(r1)
                merged = True
                break
        if not merged:
            output.append([r1])
    return output

def iter_gmap_sam_for_fusion(gmap_sam_filename, fusion_candidates, transfrag_len_dict):
    """
    Iterate through a sorted GMAP SAM file
    Continuously yield a group of overlapping records {'+': [r1, r2, ...], '-': [r3, r4....]}
    """
    records = []
    iter = BioReaders.GMAPSAMReader(gmap_sam_filename, True, query_len_dict=transfrag_len_dict)
    for r in iter:
        if r.qID in fusion_candidates:
            records = [r]
            break

    for r in iter:
        if len(records) >= 1 and (r.sID==records[-1].sID and r.sStart < records[-1].sStart):
            print("ERROR: SAM file is NOT sorted. ABORT!", file=sys.stderr)
            sys.exit(-1)
        if len(records) >= 1 and (r.sID != records[0].sID or r.sStart > records[-1].sEnd):
            yield(sep_by_strand(records))
            records = []
        if r.qID in fusion_candidates:
            records.append(r)

    if len(records) > 0:
        yield(sep_by_strand(records))


def find_fusion_candidates(sam_filename, query_len_dict, min_locus_coverage=.05, min_locus_coverage_bp=1, min_total_coverage=.99, min_dist_between_loci=10000, min_identity=0.95):
    """
    Return dict of
       fusion candidate qID --> list (in order) of the fusion ranges (ex: (chr3,100,200), (chr1,500,1000))
    (1) must map to 2 or more loci
    (2) minimum coverage for each loci is 5% AND minimum coverage in bp is >= 1 bp
    (3) total coverage is >= 95%
    (4) distance between the loci is at least 10kb
    """
    TmpRec = namedtuple('TmpRec', ['qCov', 'qLen', 'qStart', 'qEnd', 'sStart', 'sEnd', 'iden', 'chrom'])
    def total_coverage(tmprecs):
        tree = ClusterTree(0, 0)
        for r in tmprecs: tree.insert(r.qStart, r.qEnd, -1)
        return sum(reg[1]-reg[0] for reg in tree.getregions())

    d = defaultdict(lambda: [])
    reader = BioReaders.GMAPSAMReader(sam_filename, True, query_len_dict=query_len_dict)
    for r in reader:
        if r.sID == '*': continue
        if r.flag.strand == '+':
            d[r.qID].append(TmpRec(qCov=r.qCoverage, qLen=r.qLen, qStart=r.qStart, qEnd=r.qEnd, sStart=r.sStart, sEnd=r.sEnd, iden=r.identity, chrom=r.sID))
        else:
            d[r.qID].append(TmpRec(qCov=r.qCoverage, qLen=r.qLen, qStart=r.qLen-r.qEnd, qEnd=r.qLen-r.qStart, sStart=r.sStart, sEnd=r.sEnd, iden=r.identity, chrom=r.sID))
    fusion_candidates = {}
    for k, data in d.items():
        if len(data) > 1 and \
            all(a.iden>=min_identity for a in data) and \
            all(a.qCov>=min_locus_coverage for a in data) and \
            all(a.qCov*a.qLen >= min_locus_coverage_bp for a in data) and \
            total_coverage(data)*1./data[0].qLen >= min_total_coverage and \
            all(max(a.sStart,b.sStart)-min(a.sEnd,b.sEnd)>=min_dist_between_loci \
                           for a,b in itertools.combinations(data, 2)):
                    data.sort(key=lambda x: x.qStart)
                    #pdb.set_trace()
                    fusion_candidates[k] = [(a.chrom,a.sStart,a.sEnd) for a in data]
    return fusion_candidates

def fusion_main(fa_or_fq_filename, sam_filename, output_prefix, cluster_report_csv=None,
                is_fq=False, allow_extra_5_exons=True, skip_5_exon_alt=True,
                min_locus_coverage=.05, min_total_coverage=.99,
                min_locus_coverage_bp=1, min_dist_between_loci=10000,
                min_identity=0.95,
                is_flnc=False):
    """
    (1) identify fusion candidates (based on mapping, total coverage, identity, etc)
    (2) group/merge the fusion exons, using an index to point to each individual part
    (3) use BranchSimple to write out a tmp GFF where
         PBfusion.1.1 is the first part of a fusion gene
         PBfusion.1.2 is the second part of a fusion gene
    (4) read the tmp file from <3> and modify it so that
         PBfusion.1 just represents the fusion gene (a single transcript GFF format)
    """
    compressed_records_pointer_dict = defaultdict(lambda: [])
    merged_exons = []
    merged_i = 0

    # step (0). check for duplicate IDs
    check_ids_unique(fa_or_fq_filename, is_fq=is_fq)

    # step (1). identify fusion candidates
    bs = branch_simple2.BranchSimple(fa_or_fq_filename, is_fq=is_fq)
    fusion_candidates = find_fusion_candidates(sam_filename, bs.transfrag_len_dict, min_locus_coverage, min_locus_coverage_bp, min_total_coverage, min_dist_between_loci, min_identity=min_identity)

    # step (2). merge the fusion exons
    for recs in iter_gmap_sam_for_fusion(sam_filename, fusion_candidates, bs.transfrag_len_dict):
        for v in recs.values():
            if len(v) > 0:
                o = merge_fusion_exons(v, max_fusion_point_dist=100, max_exon_end_dist=0, allow_extra_5_exons=allow_extra_5_exons)
                for group in o:
                    merged_exons.append(group)
                    for r in group: compressed_records_pointer_dict[r.qID].append(merged_i)
                    merged_i += 1

    # step (3). use BranchSimple to write a temporary file
    f_group = open('branch_tmp.group.txt', 'w')
    gene_index = 1
    already_seen = set()
    for qid,indices in compressed_records_pointer_dict.items():
        combo = tuple(indices)
        if combo in already_seen:
            #print("combo seen:", combo)
            continue
        already_seen.add(combo)

        # print the fusion transcript in order of the transcription
        for i in indices:
            bs.cuff_index = gene_index # for set to the same
            records = merged_exons[i]
            isoform_index = get_isoform_index(fusion_candidates[qid], records[0].sID, records[0].sStart, records[0].sEnd)
            f_group.write("{p}.{i}.{j}\t{ids}\n".format(p="PBfusion", i=gene_index, j=isoform_index+1, ids=",".join(r.qID for r in records)))
        gene_index += 1
    f_group.close()

    # step (4). read the tmp file and modify to display per fusion gene
    # IMPORTANT: sometimes a fusion can involve more than 2 loci!
    f_group = open(output_prefix + '.group.txt', 'w')
    group_info = {} # ex: PBfusion.1 --> [id1, id2, id3...]
    count = 0
    with open('branch_tmp.group.txt') as f:
        while True:
            line = f.readline().strip()
            if len(line) == 0: break
            pbid1, groups1 = line.strip().split('\t')
            group = set(groups1.split(','))
            while True:
                cur_pos = f.tell()
                line = f.readline().strip()
                if len(line) == 0: break
                new_pbid, new_group = line.strip().split('\t')
                if new_pbid.split('.')[1]!=pbid1.split('.')[1]:
                    f.seek(cur_pos)
                    break
                else: # still in the same fusion group
                    group = group.intersection(new_group.split(','))
            f_group.write("{0}\t{1}\n".format(pbid1[:pbid1.rfind('.')], ",".join(group)))
            group_info[pbid1[:pbid1.rfind('.')]] = list(group)
            count += 1
    f_group.close()
    #os.remove('branch_tmp.group.txt')

    gff_filename = output_prefix + '.gff'
    group_filename = output_prefix + '.group.txt'
    if is_fq:
        output_filename = output_prefix + '.rep.fq'
    else:
        output_filename = output_prefix + '.rep.fa'
    pick_rep(fa_or_fq_filename, sam_filename, gff_filename, group_filename, output_filename, fusion_candidates, is_fq=is_fq)

    print("{0} fusion candidates identified.".format(count), file=sys.stdout)
    print("Output written to: {0}.gff, {0}.group.txt, {1}".format(output_prefix, output_filename), file=sys.stdout)

    # (optional) step 5. get count information
    if cluster_report_csv is not None:
        get_abundance_post_collapse(output_prefix, cluster_report_csv, output_prefix)
        print("Count information written to: {0}.abundance.txt".format(output_prefix), file=sys.stdout)
    elif is_flnc:
        print("Input is FLNC. Outputting FLNC counts per fusion.")
        with open(output_prefix+'.abundance.txt', 'w') as f:
            f.write("pbid\tcount_fl\n")
            for pbid, members in group_info.items():
                f.write("{0}\t{1}\n".format(pbid, len(members)))
        print("Count information written to: {0}.abundance.txt".format(output_prefix), file=sys.stdout)



if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("--input", help="Input FA/FQ filename")
    parser.add_argument("--fq", default=False, action="store_true", help="Input is a fastq file (default is fasta)")
    parser.add_argument("-s", "--sam", required=True, help="Sorted GMAP SAM filename")
    parser.add_argument("-o", "--prefix", required=True, help="Output filename prefix")
    parser.add_argument("--cluster_report_csv", help="cluster_report.csv, optional, if given will generate count info.")
    parser.add_argument("--is_flnc", action="store_true", default=False, help="Input are individual FLNC reads. If this option used, --cluster_report_csv should not be called. An FLNC-based count file will be output. (default: off)")
    parser.add_argument("--dun-merge-5-shorter", action="store_false", dest="allow_extra_5exon", default=True, help="Don't collapse shorter 5' transcripts (default: turned off)")
    parser.add_argument("-c", "--min_locus_coverage", type=float, default=0.05, help="Minimum per-locus coverage in percentage (default: 0.05)")
    parser.add_argument("--min_locus_coverage_bp", type=int, default=1, help="Minimum per-locus coverage in bp (default: 1 bp)")
    parser.add_argument("-t", "--min_total_coverage", type=float, default=0.99, help="Minimum total coverage (default: 0.99)")
    parser.add_argument("-d", "--min_dist_between_loci", type=int, default=10000, help="Minimum distance between loci, in bp (default: 10000)")
    parser.add_argument("-i", "--min_identity", type=float, default=0.95, help="Minimum alignment identity (default: 0.95)")

    args = parser.parse_args()

    assert 0 < args.min_identity <= 1

    if args.is_flnc and args.cluster_report_csv is not None:
        print("ERROR: Cannot use both --is_flnc and --cluster_report_csv at the same time. Abort!", file=sys.stderr)
        sys.exit(-1)

    fusion_main(args.input, args.sam, args.prefix, args.cluster_report_csv,
                is_fq=args.fq, allow_extra_5_exons=args.allow_extra_5exon,
                skip_5_exon_alt=False,
                min_locus_coverage=args.min_locus_coverage, min_locus_coverage_bp=args.min_locus_coverage_bp,
                min_total_coverage=args.min_total_coverage,
                min_dist_between_loci=args.min_dist_between_loci,
                min_identity=args.min_identity,
                is_flnc=args.is_flnc)

