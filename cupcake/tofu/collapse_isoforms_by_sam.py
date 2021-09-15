#!/usr/bin/env python
__author__ = 'etseng@pacb.com'

"""
Identical to the collapse script provided in cDNA_primer (ToFU) GitHub.

Takes a SAM file and the input fasta/fastq used to produce the SAM file,
filter out alignments based on low coverage/identity and collapse/merge
any identical isoforms based on the aligned exonic structure.

Example:
collapse_isoforms_by_sam.py --input test.fq --fq -s test.fq.sorted.sam --dun-merge-5-shorter -o test

Suggested scripts to follow up with:
   get_abundance_post_collapse.py: create count (absolute and normalized) information post collapse
   filter_by_count.py: filter away based on FL count support
   filter_away_subset.py (if collapse is run with --dun-merge-5-shorter)
"""

import os, sys
from collections import defaultdict
import threading
from multiprocessing import Process

from bx.intervals import IntervalTree
from Bio import SeqIO
from Bio.SeqIO import SeqRecord
from Bio.Seq import Seq

from cupcake.tofu.utils import check_ids_unique
from cupcake.tofu.branch import branch_simple2
from cupcake.tofu import compare_junctions
from cupcake.io import GFF

import pysam



def pick_rep(fa_fq_filename, gff_filename, group_filename, output_filename, is_fq=False, pick_least_err_instead=False, bad_gff_filename=None, fafq_dict=None):
    """
    For each group, select the representative record

    If is FASTA file (is_fa False) -- then always pick the longest one
    If is FASTQ file (is_fq True) -- then
          If pick_least_err_instead is True, pick the one w/ least number of expected base errors
          Else, pick the longest one
    """
    if fafq_dict is None:
        fd = SeqIO.to_dict(SeqIO.parse(open(fa_fq_filename), 'fastq' if is_fq else 'fasta'))
    else:
        fd = fafq_dict
    fout = open(output_filename, 'w')

    coords = {}
    for r in collapseGFFReader(gff_filename):
        tid = r.seqid
        coords[tid] = "{0}:{1}-{2}({3})".format(r.chr, r.start, r.end, r.strand)

    if bad_gff_filename is not None:
        for r in collapseGFFReader(bad_gff_filename):
            tid = r.seqid
            coords[tid] = "{0}:{1}-{2}({3})".format(r.chr, r.start, r.end, r.strand)
            
    for line in open(group_filename):
        pb_id, members = line.strip().split('\t')
        print("Picking representative sequence for {0}".format(pb_id), file=sys.stdout)
        best_rec = None
        #best_id = None
        #best_seq = None
        #best_qual = None
        best_err = 9999999
        err = 9999999
        max_len = 0
        for x in members.split(','):
            if is_fq and pick_least_err_instead:
                err = sum(i**-(i/10.) for i in fd[x].letter_annotations['phred_quality'])
            if (is_fq and pick_least_err_instead and err < best_err) or ((not is_fq or not pick_least_err_instead) and len(fd[x].seq) >= max_len):
                best_rec = fd[x]
                #best_id = x
                #best_seq = fd[x].seq
                #if is_fq:
                #    best_qual = fd[x].quality
                #    best_err = err
                max_len = len(fd[x].seq)

        _id_ = "{0}|{1}|{2}".format(pb_id, coords[pb_id], best_rec.id)
        best_rec.id = _id_
        SeqIO.write(best_rec, fout, 'fastq' if is_fq else 'fasta')

    fout.close()

def collapse_fuzzy_junctions(gff_filename, group_filename, allow_extra_5exon, internal_fuzzy_max_dist):
    def can_merge(m, r1, r2):
        if m == 'exact':
            return True
        else:
            if not allow_extra_5exon:
                return False
        # below is continued only if (a) is 'subset' or 'super' AND (b) allow_extra_5exon is True
        if m == 'subset':
            r1, r2 = r2, r1 #  rotate so r1 is always the longer one
        if m == 'super' or m == 'subset':
            n2 = len(r2.ref_exons)
            # check that (a) r1 and r2 end on same 3' exon, that is the last acceptor site agrees
            # AND (b) the 5' start of r2 is sandwiched between the matching r1 exon coordinates
            if r1.strand == '+':
                return abs(r1.ref_exons[-1].start - r2.ref_exons[-1].start) <= internal_fuzzy_max_dist and \
                    r1.ref_exons[-n2].start <= r2.ref_exons[0].start < r1.ref_exons[-n2].end
            else:
                return abs(r1.ref_exons[0].end - r2.ref_exons[0].end) <= internal_fuzzy_max_dist and \
                    r1.ref_exons[n2-1].start <= r2.ref_exons[-1].end < r1.ref_exons[n2-1].end
        return False

    d = {}
    recs = defaultdict(lambda: {'+':IntervalTree(), '-':IntervalTree()}) # chr --> strand --> tree
    fuzzy_match = defaultdict(lambda: [])
    for r in GFF.collapseGFFReader(gff_filename):
        d[r.seqid] = r
        has_match = False
        r.segments = r.ref_exons
        for r2 in recs[r.chr][r.strand].find(r.start, r.end):
            r2.segments = r2.ref_exons
            m = compare_junctions.compare_junctions(r, r2, internal_fuzzy_max_dist=internal_fuzzy_max_dist, max_5_diff=args.max_5_diff, max_3_diff=args.max_3_diff)
            if can_merge(m, r, r2):
                fuzzy_match[r2.seqid].append(r.seqid)
                has_match = True
                break
        if not has_match:
            recs[r.chr][r.strand].insert(r.start, r.end, r)
            fuzzy_match[r.seqid] = [r.seqid]

    group_info = {}
    with open(group_filename) as f:
        for line in f:
            pbid, members = line.strip().split('\t')
            group_info[pbid] = [x for x in members.split(',')]

    # pick for each fuzzy group the one that has the most exons
    keys = list(fuzzy_match.keys())
    keys.sort(key=lambda x: int(x.split('.')[1]))
    f_gff = open(gff_filename+'.fuzzy', 'w')
    f_group = open(group_filename+'.fuzzy', 'w')
    for k in keys:
        all_members = []
        best_pbid, best_size, best_num_exons = fuzzy_match[k][0], len(group_info[fuzzy_match[k][0]]), len(d[fuzzy_match[k][0]].ref_exons)
        all_members += group_info[fuzzy_match[k][0]]
        for pbid in fuzzy_match[k][1:]:
            _num_exons = len(d[pbid].ref_exons)
            _size = len(group_info[pbid])
            all_members += group_info[pbid]
            if _num_exons > best_num_exons or (_num_exons==best_num_exons and _size>best_size):
                best_pbid, best_size, best_num_exons = pbid, _size, _num_exons
        GFF.write_collapseGFF_format(f_gff, d[best_pbid])
        f_group.write("{0}\t{1}\n".format(best_pbid, ",".join(all_members)))
    f_gff.close()
    f_group.close()

    return fuzzy_match

def multiprocess_predefine_regions(aligned_bam_filename, n_chunks):
    """
    Read through an aligned sorted BAM file and split it into regions
    grouping them properly by overlaps (in this version we are strand-insensitive)
    :return: list of [start_index, end_index) of non-overolapping regions,
             CHUNKED list of [start_index, end_index) based on <n> chunks given
    """
    reader = pysam.AlignmentFile(open(aligned_bam_filename), 'rb', check_sq=False)
    prev = None
    for prev in reader:
        if not prev.is_unmapped:
            break
    if prev is None:
        print(f"BAM file {aligned_bam_filename} does not contain any mapped records. ABORT!", file=sys.stderr)
        sys.exit(-1)
    max_end = prev.reference_end
    prev_index = 0
    cur_index = 0
    region_list = [] # list of regions defined by [start_index, end_index)
    for r in reader:
        cur_index += 1
        if r.is_unmapped: continue
        if r.reference_name == prev.reference_name and r.reference_start < prev.reference_start:
            print(f"BAM file {aligned_bam_filename }is NOT sorted. ABORT!", file=sys.stderr)
            sys.exit(-1)
        if r.reference_name != prev.reference_name or r.reference_start > max_end:
            # we've found a new non-overlapping region
            region_list.append((prev_index, cur_index))
            prev = r
            prev_index = cur_index
            max_end = prev.reference_end
        else:
            max_end = max(max_end, r.reference_end)
    if cur_index > prev_index:
        region_list.append((prev_index, cur_index+1))

    # split the regions by cpu given
    chunk_size = len(region_list) // n_chunks
    chunked_region_index_list = []
    for i in range(n_chunks-1):
        start_index = region_list[chunk_size * i][0]  # the (0-based) start index of the first region
        end_index = region_list[chunk_size * (i + 1)][0]  # the (1-based) end index of the last region
        chunked_region_index_list.append((start_index, end_index))
    chunked_region_index_list.append((region_list[chunk_size*(n_chunks-1)][0],
                                      region_list[-1][1]))

    return region_list, chunked_region_index_list

def multiprocess_helper(start_index, end_index, args, cov_threshold, f_good, f_bad, f_txt, ignored_fout):
    b = branch_simple2.BranchSimple(args.input,
                                    cov_threshold=cov_threshold,
                                    min_aln_coverage=args.min_aln_coverage,
                                    min_aln_identity=args.min_aln_identity,
                                    is_fq=args.fq,
                                    max_5_diff=args.max_5_diff,
                                    max_3_diff=args.max_3_diff)

    iter = b.iter_gmap_sam(args.bam, ignored_fout, type='BAM', bam_start_index=start_index, bam_end_index=end_index)
    for recs in iter: # recs is {'+': list of list of records, '-': list of list of records}
        for v in recs.values():
            for v2 in v:
                if len(v2) > 0: b.process_records(v2, args.allow_extra_5exon, False, f_good, f_bad, f_txt)

    # for multiprocessing, we close the files here!
    f_good.close()
    if f_bad!=f_good: f_bad.close()
    f_txt.close()
    ignored_fout.close()


import re
from cupcake.io.GFF import collapseGFFReader, write_collapseGFF_format
pbid_rex = re.compile('PB.(\d+).(\d+)')

def multiprocess_combine_result(in_gff_list, f_gff, in_group_list, f_group):
    # first GFF we just directly output it, no offset needed
    reader_gff = collapseGFFReader(in_gff_list[0])
    reader_group = open(in_group_list[0]) if in_group_list is not None else None
    for r in reader_gff:
        write_collapseGFF_format(f_gff, r)
        if reader_group is not None:
            last_pbid, members = reader_group.readline().strip().split()
            assert last_pbid == r.seqid
            f_group.write(f"{last_pbid}\t{members}\n")


    pbid_offset = int(pbid_rex.match(last_pbid).group(1))

    if in_group_list is None: in_group_list = [None]*len(in_gff_list)

    for (in_gff, in_group) in zip(in_gff_list[1:], in_group_list[1:]):
        print("Combining {0} into {1} with offset PB.{2}...".format(in_gff, f_gff.name, pbid_offset))
        reader_gff = collapseGFFReader(in_gff)
        reader_group = open(in_group, 'r')
        for r in reader_gff:
            m = pbid_rex.match(r.seqid)
            gene_index, isoform_index = int(m.group(1)), m.group(2)
            last_pbid = gene_index + pbid_offset
            new_geneid = "PB.{0}".format(last_pbid)
            new_pbid = "PB.{0}.{1}".format(last_pbid, isoform_index)
            if reader_group is not None:
                pbid, members = reader_group.readline().strip().split()
                assert pbid == r.seqid
                f_group.write(f"{new_pbid}\t{members}\n")
            # write the GFF record, updating the PBID
            r.seqid = new_pbid
            r.geneid = new_geneid
            write_collapseGFF_format(f_gff, r)
        pbid_offset = last_pbid


def main(args):

    ### sanity check that input file and input SAM exists
    if args.input is None:
        if args.bam is None:
            print("ERROR! --bam must be provided if input fasta is not given. Abort!")
            sys.exit(-1)
        print("WARNING: input fasta not given. Will use query sequence in the BAM file.")
    else:
        if not os.path.exists(args.input):
            print("Input file {0} does not exist. Abort.".format(args.input), file=sys.stderr)
            sys.exit(-1)

    if args.sam is not None and not os.path.exists(args.sam):
        print("SAM file {0} does not exist. Abort.".format(args.sam), file=sys.stderr)
        sys.exit(-1)

    if args.bam is not None and not os.path.exists(args.bam):
        print("BAM file {0} does not exist. Abort.".format(args.bam), file=sys.stderr)
        sys.exit(-1)

    # check for duplicate IDs
    if args.input is not None:
        check_ids_unique(args.input, is_fq=args.fq)

    ignored_fout = open(args.prefix + '.ignored_ids.txt', 'w')

    if args.flnc_coverage > 0:
        f_good = open(args.prefix + '.collapsed.good.gff', 'w')
        f_bad = open(args.prefix + '.collapsed.bad.gff', 'w')
        cov_threshold = args.flnc_coverage
    else:
        f_good = open(args.prefix + '.collapsed.gff', 'w')
        f_bad = f_good
        cov_threshold = 1
    f_txt = open(args.prefix + '.collapsed.group.txt', 'w')

    if args.cpus == 1:
        b = branch_simple2.BranchSimple(args.input, cov_threshold=cov_threshold, min_aln_coverage=args.min_aln_coverage,
                                        min_aln_identity=args.min_aln_identity, is_fq=args.fq,
                                        max_5_diff=args.max_5_diff, max_3_diff=args.max_3_diff)
        if args.bam is not None:
            iter = b.iter_gmap_sam(args.bam, ignored_fout, type='BAM')
        else:
            iter = b.iter_gmap_sam(args.sam, ignored_fout, type='SAM')
        for recs in iter: # recs is {'+': list of list of records, '-': list of list of records}
            for v in recs.values():
                for v2 in v:
                    if len(v2) > 0: b.process_records(v2, args.allow_extra_5exon, False, f_good, f_bad, f_txt)
    else:
        # need to first predefine the regions
        region_list_ignore, chunk_regions_list = multiprocess_predefine_regions(args.bam, args.cpus)
        assert len(chunk_regions_list) == args.cpus

        if args.flnc_coverage > 0:
            f_good_pool = [open(args.prefix + '.collapsed.good.gff' + str(i), 'w') for i in range(args.cpus)]
            f_bad_pool = [open(args.prefix + '.collapsed.bad.gff' + str(i), 'w') for i in range(args.cpus)]
        else:
            f_good_pool = [open(args.prefix + '.collapsed.gff' + str(i), 'w') for i in range(args.cpus)]
            f_bad_pool = f_good_pool
        f_txt_pool = [open(args.prefix + '.collapsed.group.txt' + str(i), 'w') for i in range(args.cpus)]
        f_ignore_pool = [open(args.prefix + '.ignords_ids.txt' + str(i), 'w') for i in range(args.cpus)]

        pool = []
        for i in range(args.cpus):
            p = Process(target=multiprocess_helper,
                                 args=(chunk_regions_list[i][0], chunk_regions_list[i][1], args, cov_threshold,
                                       f_good_pool[i], f_bad_pool[i], f_txt_pool[i], f_ignore_pool[i], ))
            p.start()
            pool.append(p)
            # NOTE: f_good_pool[i] and f_bad_pool[i] and f_txt_pool[i] actually will get file CLOSED
        for p in pool:
            p.join()

        # for .ignore_ids.txt we can just concatenate
        for f in f_ignore_pool:
            with open(f.name, 'r') as h:
                ignored_fout.write(h.read())
        multiprocess_combine_result([f.name for f in f_good_pool], f_good, [f.name for f in f_txt_pool], f_txt)
        if args.flnc_coverage > 0:
            multiprocess_combine_result([f.name for f in f_bad_pool], f_bad, None, None)

        # now we can delete the chunked results
        for f in f_good_pool:
            os.remove(f.name)
        if f_good_pool!=f_bad_pool:
            for f in f_bad_pool:
                os.remove(f.name)
        for f in f_ignore_pool:
            os.remove(f.name)
        for f in f_txt_pool:
            os.remove(f.name)


    ignored_fout.close()
    f_good.close()
    f_bad.close()
    f_txt.close()

    if args.max_fuzzy_junction > 0: # need to further collapse those that have fuzzy junctions!
        collapse_fuzzy_junctions(f_good.name, f_txt.name, args.allow_extra_5exon, internal_fuzzy_max_dist=args.max_fuzzy_junction)
        os.rename(f_good.name, f_good.name+'.unfuzzy')
        os.rename(f_txt.name, f_txt.name+'.unfuzzy')
        os.rename(f_good.name+'.fuzzy', f_good.name)
        os.rename(f_txt.name+'.fuzzy', f_txt.name)

    if args.fq:
        outfile = args.prefix+".collapsed.rep.fq"
    else:
        outfile = args.prefix+".collapsed.rep.fa"

    if args.input is not None:
        is_fq = args.fq
        fafq_dict = None
    else:
        is_fq = False
        fafq_dict = {}
        for r in pysam.AlignmentFile(open(args.bam), 'rb', check_sq=False):
            fafq_dict[r.qname] = SeqRecord(Seq(r.query), id=r.qname)

    if args.allow_extra_5exon: # 5merge, pick longest
        pick_rep(args.input, f_good.name, f_txt.name, outfile,
                 is_fq=is_fq,
                 pick_least_err_instead=False,
                 bad_gff_filename=f_bad.name,
                 fafq_dict=fafq_dict)
    else:
        pick_rep(args.input, f_good.name, f_txt.name, outfile,
                 is_fq=is_fq,
                 pick_least_err_instead=True,
                 bad_gff_filename=f_bad.name,
                 fafq_dict=fafq_dict)

    if args.gen_mol_count:
        outfile = args.prefix + '.collapsed.abundance.txt'
        with open(outfile, 'w') as f:
            f.write("pbid\tcount_fl\n")
            for line in open(f_txt.name):
                pbid, members = line.strip().split()
                f.write("{0}\t{1}\n".format(pbid, members.count(',')+1))


    print("Ignored IDs written to: {0}".format(ignored_fout.name), file=sys.stdout)
    print("Output written to: {0}\n{1}\n{2}\n{3}\n".format(f_good.name, f_txt.name, outfile, args), file=sys.stdout)

if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("--input", help="Input FA/FQ filename")
    parser.add_argument("--fq", default=False, action="store_true", help="Input is a fastq file (default is fasta)")
    parser.add_argument("-s", "--sam", required=False, default=None, help="Sorted SAM filename, use either this or --bam")
    parser.add_argument("-b", "--bam", required=False, default=None, help="Sorted BAM filename, use either this or --sam")
    parser.add_argument("-o", "--prefix", required=True, help="Output filename prefix")
    parser.add_argument("-c", "--min-coverage", dest="min_aln_coverage", type=float, default=.99, help="Minimum alignment coverage (default: 0.99)")
    parser.add_argument("-i", "--min-identity", dest="min_aln_identity", type=float, default=.95, help="Minimum alignment identity (default: 0.95)")
    parser.add_argument("--max_fuzzy_junction", default=5, type=int, help="Max fuzzy junction dist (default: 5 bp)")
    parser.add_argument("--max_5_diff", default=1000, type=int, help="Maximum allowed 5' difference if on same exon (default: 1000 bp)")
    parser.add_argument("--max_3_diff", default=100, type=int, help="Maximum allowed 3' difference if on same exon (default: 100 bp)")
    parser.add_argument("--flnc_coverage", dest="flnc_coverage", type=int, default=-1, help="Minimum # of FLNC reads, only use this for aligned FLNC reads, otherwise results undefined!")
    parser.add_argument("--gen_mol_count", action="store_true", default=False, help="Generate a .abundance.txt file based on the number of input sequences collapsed. Use only if input is FLNC or UMI-dedup output (default: off)")
    parser.add_argument("--dun-merge-5-shorter", action="store_false", dest="allow_extra_5exon", default=True, help="Don't collapse shorter 5' transcripts (default: turned off)")
    parser.add_argument("--cpus", default=1, type=int, help="Number of CPUs for parallelization (default: 1)")

    args = parser.parse_args()

    if args.sam and args.bam:
        print("Must provide only a SAM via --sam or only a BAM via --bam but not both! Abort!")
        sys.exit(-1)
    if args.sam is None and args.bam is None:
        print("Must provide an input aligned SAM via --sam or aligned BAM via --bam. Abort!")
        sys.exit(-1)

    main(args)


