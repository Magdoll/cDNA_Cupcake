#!/usr/bin/env python
__author__ = 'etseng@pacb.com'

"""
Should be used *after* summarize_sample_GFF_junctions.py is called to generate junction report.

Looks through the junction reports and scrub it, retaining only junctions that
 meet one or more of the following criteria:

(a). be the only junction in that "label" <-- i.e. you must keep this junction because there's no other similar ones to it
(b). have previous annotation (annotation='Y') <-- always trust reference annotation
(c). have at least X sample supporting it (`num_sample>=T` or `num_transcript>=S`) where `S`, `T` is user-defined
(d). additionally, --accept-all-canonical <-- meaning all canonical junctions are accepted

Output: scrubbed.junctions.bed
"""
import os, sys
from collections import defaultdict
from csv import DictReader, DictWriter
from bx.intervals import IntervalTree, Interval
from Bio import SeqIO

import cupcake.io.GFF as GFF
import cupcake.tofu.counting.chain_samples as sp

fields_to_add = ['count_fl','count_nfl','count_nfl_amb','norm_fl','norm_nfl','norm_nfl_amb']

def read_junction_report(filename):
    """
    tab-delimited with header:
           chr     left    right   strand  num_transcript  num_sample      genome  annotation      label

    return: dict of label --> records with that label
    """
    reader = DictReader(open(filename), delimiter='\t')
    r = next(reader)

    cur_label, cur = r['label'], [r]
    for r in reader:
        if r['label']!=cur_label:
            yield cur_label, cur
            cur_label, cur = r['label'], [r]
        else:
            cur.append(r)
    yield cur_label, cur

def scrub_junction_by_label(junctions, min_sample=2, min_transcript=2, accept_all_canonical=False):
    """
    input: a list of junctions (records) that are of the same label
    output: a list of "scrubbed" junctions
    """
    def pass_criteria(r):
        if int(r['num_sample']) >= min_sample: return True
        if int(r['num_transcript']) >= min_transcript: return True
        if accept_all_canonical and r['genome']=='GT-AG': return True
        if r['annotation']=='Y': return True

    if len(junctions)==1:
        # only one junction in this label, must accept!
        return junctions

    good = []
    for r in junctions:
        if pass_criteria(r):
            good.append(r)

    if len(good)==0: # nothing passed criteria! must pick one. we will pick one that has the highest transcript support
        junctions.sort(key=lambda r: int(r['num_transcript']), reverse=True)
        good.append(junctions[0])
    return good


def find_best_match_junction(tree, donor, accep, max_diff=20):
    """
    donor, accept -- both should be 0-based
    """
    hits = tree.find(donor, accep)
    if len(hits) == 0:
        return None
    elif len(hits) == 1:
        if hits[0].start-donor > max_diff or hits[0].end-accep > max_diff:
            return None
        return hits[0]
    else: # multiple hits, find the closest one
        diff = []
        for h in hits:
            if h.start-donor > max_diff or h.end-accep > max_diff: continue
            diff.append((abs(h.start-donor)+abs(h.end-accep), h))
        diff.sort(key=lambda x: x[0])
        return diff[0][1]

def scrub_ref_exons(r, tree):
    n = len(r.ref_exons)
    new_ref_exons = []
    cur_start = r.ref_exons[0].start
    for i in range(n-1):
        donor = r.ref_exons[i].end-1 # make it 0-based
        accep = r.ref_exons[i+1].start # start is already 0-based
        match = find_best_match_junction(tree[r.chr, r.strand], donor, accep)
        if match is None:
            print("donor-acceptor site {0},{1},{2}-{3} has no hit in tree!".format(\
                r.chr, r.strand, donor, accep), file=sys.stderr)
            return None

        new_ref_exons.append(Interval(cur_start, match.start+1))
        cur_start = match.end
    new_ref_exons.append(Interval(cur_start, r.ref_exons[-1].end))
    return new_ref_exons

def read_scrubbed_junction_to_tree(junction_filename):
    tree = defaultdict(lambda: IntervalTree())
    f = open(junction_filename)
    if not f.readline().startswith('track'): f.seek(0)
    for line in f:
        raw = line.strip().split('\t')
        if len(raw) == 4: chrom, left, right, strand = raw
        elif len(raw) == 6: chrom, left, right, _name, _count, strand = raw
        else: raise Exception("Expects junction BED file to have either 4 or 6 columns! Saw {0}!".format(len(raw)))
        left, right = int(left), int(right) # already 0-based start, 0-based end
        tree[chrom,strand].add(left, right, Interval(left, right))
    return tree

def scrub_junctions(report_filename, output_filename, min_sample, min_transcript, accept_all_canonical):
    tree = defaultdict(lambda: IntervalTree())
    f = open(output_filename, 'w')
    for _label, junctions in read_junction_report(report_filename):
        good = scrub_junction_by_label(junctions, min_sample, min_transcript, accept_all_canonical)
        for r in good:
            a, b = int(r['left']), int(r['right']) # 0-based start, 0-basde end
            f.write("{chrom}\t{left}\t{right}\t{strand}\n".format(\
                chrom=r['chr'], left=r['left'], right=r['right'], strand=r['strand']))
            tree[r['chr'],r['strand']].add(a, b, Interval(a, b))
    f.close()
    return tree

def scrub_sample_GFFs(sample_dirs, gff_filename, count_filename, group_filename, fastq_filename, output_prefix, tree):

    for sample_name, d in sample_dirs.items():
        outf = open(os.path.join(d, output_prefix+'.gff.tmp'), 'w')
        for r in GFF.collapseGFFReader(os.path.join(d, gff_filename)):
            n = len(r.ref_exons)
            if n == 1:
                GFF.write_collapseGFF_format(outf, r)

            new_ref_exons = scrub_ref_exons(r, tree)
            if new_ref_exons is None:
                print("No changes made due to error:", r.seqid, file=sys.stderr)
            else:
                #print "before:", r.ref_exons
                #print "after :", new_ref_exons
                r.ref_exons = new_ref_exons
            GFF.write_collapseGFF_format(outf, r)
        outf.close()
        cleanup_scrubbed_files_redundancy(outf.name, \
                                          os.path.join(d, group_filename), \
                                          os.path.join(d, count_filename), \
                                          os.path.join(d, fastq_filename) if fastq_filename is not None else None,
                                          os.path.join(d, output_prefix))


def read_count_file(count_filename):
    f = open(count_filename)
    count_header = ''
    while True:
        cur_pos = f.tell()
        line = f.readline()
        if not line.startswith('#'):
            f.seek(cur_pos)
            break
        else:
            count_header += line
    d = dict((r['pbid'], r) for r in DictReader(f, delimiter='\t'))
    f.close()
    return d, count_header

def read_group_file(group_filename):
    group_info = {} # key: PB.1.1 --> list of HQ isoforms
    for line in open(group_filename):
        a, b = line.strip().split('\t')
        group_info[a] = b
    return group_info

def cleanup_scrubbed_files_redundancy(gff_filename, group_filename, count_filename, fastq_filename, output_prefix):

    junction_seen = defaultdict(lambda: defaultdict(lambda: [])) # key (chr,strand) -> dict of (series of junctions) -> record
    for r in GFF.collapseGFFReader(gff_filename):
        n = len(r.ref_exons)
        if n == 1:
            junc_str = str(r.start)+','+str(r.end)
            junction_seen[r.chr, r.strand][junc_str] = [r]
        else:
            junc_str = ",".join(str(r.ref_exons[i].end)+','+str(r.ref_exons[i+1].start) for i in range(n-1))
            junction_seen[r.chr, r.strand][junc_str].append(r)

    # write out cleaned GFF
    outf = open(output_prefix+'.gff', 'w')
    outf2 = open(output_prefix+'.merged_ids.txt', 'w')
    merged = {}
    keys = list(junction_seen.keys())
    keys.sort()
    for k in keys:
        for bunch in junction_seen[k].values():
            if len(bunch) == 1: # just one record, write it out
                r = bunch[0]
                GFF.write_collapseGFF_format(outf, r)
                merged[r.seqid] = [r.seqid]
            else:
                # find the representative
                r = bunch[0]
                for r2 in bunch[1:]:
                    if r2.end-r2.start > r.end-r.start:
                        r = r2
                GFF.write_collapseGFF_format(outf, r)
                merged[r.seqid] = [x.seqid for x in bunch]
            outf2.write("{0}\t{1}\n".format(r.seqid, ",".join(merged[r.seqid])))
    outf.close()
    outf2.close()

    count_d, count_header = read_count_file(count_filename)
    # write out count file
    outf = open(output_prefix + '.abundance.txt', 'w')
    outf.write(count_header)
    writer = DictWriter(outf, fieldnames=['pbid','count_fl','count_nfl','count_nfl_amb','norm_fl','norm_nfl','norm_nfl_amb'], \
                        delimiter='\t', lineterminator='\n')
    writer.writeheader()
    for pbid, bunch in merged.items():
        # combine the counts
        r = count_d[bunch[0]]
        r['pbid'] = pbid
        for field in fields_to_add:
            r[field] = float(r[field])
        for _id in bunch[1:]:
            for field in fields_to_add:
                r[field] += float(count_d[_id][field])
        writer.writerow(r)
    outf.close()

    group_info = read_group_file(group_filename)
    # write out group file
    outf = open(output_prefix + '.group.txt', 'w')
    for pbid, bunch in merged.items():
        # combine the groups
        g = [group_info[bunch[0]]]
        for _id in bunch[1:]:
            g.append(group_info[_id])
        outf.write("{0}\t{1}\n".format(pbid, ",".join(g)))
    outf.close()

    # write out fastq file if present
    if fastq_filename is not None:
        outf = open(output_prefix+'.rep.fq', 'w')
        for r in SeqIO.parse(open(fastq_filename), 'fastq'):
            if r.id.split('|')[0] in merged or r.id in merged:
                SeqIO.write(r, outf, 'fastq')
        outf.close()

    print("scrubbed files written: {0}.gff, {0}.group.txt, {0}.abundance.txt, {0}.merged_ids.txt".format(output_prefix), file=sys.stderr)



if __name__=="__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("sample_config")
    parser.add_argument("summary_report")
    parser.add_argument("output_prefix")
    parser.add_argument("-S", "--min_sample", type=int, default=1, help="Minimum number of samples as evidence (default: 1)")
    parser.add_argument("-T", "--min_transcript", type=int, default=2, help="Minimum number of transcripts as evidence (default: 2)")
    #parser.add_argument("-C", "--accept_all_canonical", action="store_true", default=False, help="Accept all canonical jucntions (default: false)")
    parser.add_argument("--scrubbed_junction_file", help="Scrubbed junction bed --- if given, directly use it to scrub GFFs.")

    args = parser.parse_args()

    sample_dirs, sample_names, group_filename, gff_filename, count_filename, fastq_filename = sp.read_config(args.sample_config)

    report_filename = args.summary_report

    if args.scrubbed_junction_file is None:
        output_filename = args.output_prefix  + ".scrubbed.junction.bed"
        tree = scrub_junctions(report_filename, output_filename, args.min_sample, args.min_transcript, True)
        print("Scrubbed junction written to: ", output_filename, file=sys.stderr)
    else:
        output_filename = args.scrubbed_junction_file
        print("Reading scrubbed junction file: ", output_filename, file=sys.stderr)
        tree = read_scrubbed_junction_to_tree(output_filename)

    scrub_sample_GFFs(sample_dirs, gff_filename, count_filename, group_filename, fastq_filename, args.output_prefix, tree)
