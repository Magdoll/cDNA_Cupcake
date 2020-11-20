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

from bx.intervals import IntervalTree
from Bio import SeqIO

from cupcake.tofu.utils import check_ids_unique
from cupcake.tofu.branch import branch_simple2
from cupcake.tofu import compare_junctions
from cupcake.io import GFF



def pick_rep(fa_fq_filename, gff_filename, group_filename, output_filename, is_fq=False, pick_least_err_instead=False, bad_gff_filename=None):
    """
    For each group, select the representative record

    If is FASTA file (is_fa False) -- then always pick the longest one
    If is FASTQ file (is_fq True) -- then
          If pick_least_err_instead is True, pick the one w/ least number of expected base errors
          Else, pick the longest one
    """
    fd = SeqIO.to_dict(SeqIO.parse(open(fa_fq_filename), 'fastq' if is_fq else 'fasta'))
    fout = open(output_filename, 'w')

    coords = {}
    for line in open(gff_filename):
        # ex: chr1    PacBio  transcript      27567   29336   .       -       .       gene_id "PB.1"; transcript_id "PB.1.1";
        raw = line.strip().split('\t')
        if raw[2] == 'transcript':
            tid = raw[-1].split('; ')[1].split()[1][1:-2]
            coords[tid] = "{0}:{1}-{2}({3})".format(raw[0], raw[3], raw[4], raw[6])

    if bad_gff_filename is not None:
        for line in open(bad_gff_filename):
            raw = line.strip().split('\t')
            if raw[2] == 'transcript':
                tid = raw[-1].split('; ')[1].split()[1][1:-2]
                coords[tid] = "{0}:{1}-{2}({3})".format(raw[0], raw[3], raw[4], raw[6])

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
                    r1.ref_exons[n2-1].start <= r2.ref_exons[-1].end < r1.ref_exons[n2].end
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


def main(args):

    ### sanity check that input file and input SAM exists
    if not os.path.exists(args.input):
        print("Input file {0} does not exist. Abort.".format(args.input), file=sys.stderr)
        sys.exit(-1)

    if not os.path.exists(args.sam):
        print("SAM file {0} does not exist. Abort.".format(args.sam), file=sys.stderr)
        sys.exit(-1)

    # check for duplicate IDs
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

    b = branch_simple2.BranchSimple(args.input, cov_threshold=cov_threshold, min_aln_coverage=args.min_aln_coverage, min_aln_identity=args.min_aln_identity, is_fq=args.fq, max_5_diff=args.max_5_diff, max_3_diff=args.max_3_diff)
    iter = b.iter_gmap_sam(args.sam, ignored_fout)
    for recs in iter: # recs is {'+': list of list of records, '-': list of list of records}
        for v in recs.values():
            for v2 in v:
                if len(v2) > 0: b.process_records(v2, args.allow_extra_5exon, False, f_good, f_bad, f_txt)

    ignored_fout.close()
    f_good.close()
    try:
        f_bad.close()
    except:
        pass
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
    if args.allow_extra_5exon: # 5merge, pick longest
        pick_rep(args.input, f_good.name, f_txt.name, outfile, is_fq=args.fq, pick_least_err_instead=False, bad_gff_filename=f_bad.name)
    else:
        pick_rep(args.input, f_good.name, f_txt.name, outfile, is_fq=args.fq, pick_least_err_instead=True, bad_gff_filename=f_bad.name)

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
    parser.add_argument("-s", "--sam", required=True, help="Sorted GMAP SAM filename")
    parser.add_argument("-o", "--prefix", required=True, help="Output filename prefix")
    parser.add_argument("-c", "--min-coverage", dest="min_aln_coverage", type=float, default=.99, help="Minimum alignment coverage (default: 0.99)")
    parser.add_argument("-i", "--min-identity", dest="min_aln_identity", type=float, default=.95, help="Minimum alignment identity (default: 0.95)")
    parser.add_argument("--max_fuzzy_junction", default=5, type=int, help="Max fuzzy junction dist (default: 5 bp)")
    parser.add_argument("--max_5_diff", default=1000, type=int, help="Maximum allowed 5' difference if on same exon (default: 1000 bp)")
    parser.add_argument("--max_3_diff", default=100, type=int, help="Maximum allowed 3' difference if on same exon (default: 100 bp)")
    parser.add_argument("--flnc_coverage", dest="flnc_coverage", type=int, default=-1, help="Minimum # of FLNC reads, only use this for aligned FLNC reads, otherwise results undefined!")
    parser.add_argument("--gen_mol_count", action="store_true", default=False, help="Generate a .abundance.txt file based on the number of input sequences collapsed. Use only if input is FLNC or UMI-dedup output (default: off)")
    parser.add_argument("--dun-merge-5-shorter", action="store_false", dest="allow_extra_5exon", default=True, help="Don't collapse shorter 5' transcripts (default: turned off)")

    args = parser.parse_args()

    main(args)


