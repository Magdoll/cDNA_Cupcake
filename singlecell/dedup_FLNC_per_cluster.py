#!/usr/bin/env python
__author__="etseng@pacb.com"

"""
After running UMI_BC_error_correct.py with --bc_rank_file and --only_top_ranked parameters,
 take the (1) .annotated.correct.csv file, and
          (2) short read cluster (cell type) csv file
          (3) fasta, gff, [faa]
          
 to de-duplicate the FLNC reads by cluster (cell type).
 
Outputs:
(1) a de-dup table of

(UMI-ed, BC-ed, gene) — pbid — associate gene — associated transcript — category — length — cluster #

(2)  A “master” file, one sequence for each [pbid] that appeared at least once in (1) 
- fasta
- gff

(3) A “per cluster” file, one sequence for each [pbid] that appeared once in each cluster
- fasta
- gff
"""

import os, sys, re
from csv import DictReader, DictWriter
from collections import Counter, defaultdict
from cupcake.io.GFF import write_collapseGFF_format, collapseGFFReader
from Bio import SeqIO

CORRECTED_CSV_FILELDS = ['id', 'pbid', 'length', 'transcript', 'gene', 'category', 'ORFgroup', 'UMI_ed', 'BC_ed']
#PBID_FORMAT = re.compile("PB.(\d+).(\d+)")

def main(corrected_csv, cluster_info, output_prefix, fasta_file=None, gff_file=None, faa_file=None):

    # read corrected CSV
    reader = DictReader(open(corrected_csv), delimiter='\t')
    for k in CORRECTED_CSV_FILELDS:
        if k not in reader.fieldnames:
            print("The following fields must exist in {0}!\n{1}".format(corrected_csv, "\n".join(CORRECTED_CSV_FILELDS)))
            sys.exit(-1)

    per_unique = {}  # tag -> record
    per_unique_count = Counter()  # tag -> number of duplicates
    per_pbid = defaultdict(lambda: {'gene':None, 'transcript':None, 'clusters':[]})  # pbid --> list of clusters it is in
    for r in reader:
        tag = "{bc}-{umi}-{gene}".format(bc=r['BC_ed'], umi=r['UMI_ed'], gene=r['gene'])
        per_unique[tag] = r
        per_unique_count[tag] += 1

    # now link barcode to cell type, also PCR dup counts
    for tag in per_unique:
        c = cluster_info[per_unique[tag]['BC_ed']]
        rec = per_unique[tag]
        rec['cluster'] = c
        rec['num_dups'] = per_unique_count[tag]
        pbid = rec['pbid']
        if pbid in per_pbid: per_pbid[pbid]['clusters'].add(c)
        else: per_pbid[pbid] = {'gene':rec['gene'], 'transcript':rec['transcript'], 'clusters':set([c])}

    # write out de-dup CSV file
    with open(output_prefix + '.csv', 'w') as f:
        writer = DictWriter(f, CORRECTED_CSV_FILELDS+['cluster', 'num_dups'], delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        keys = per_unique.keys()
        for k in sorted(keys):
            writer.writerow(per_unique[k])

    if fasta_file is not None:
        f_d = {}  # cluster --> file handle
        # writer pbid master file
        with open(output_prefix + '.fasta', 'w') as f:
            for r in SeqIO.parse(open(fasta_file), 'fasta'):
                if r.id in per_pbid:
                    newid = "{pbid}|{gene}|{transcript}|{clusters}".format(\
                            pbid=r.id,
                            gene=per_pbid[r.id]['gene'],
                            transcript=per_pbid[r.id]['transcript'],
                            clusters=";".join(per_pbid[r.id]['clusters']))
                    f.write(">{0}\n{1}\n".format(newid, r.seq))
                    for c in per_pbid[r.id]['clusters']:
                        if c not in f_d:
                            f_d[c] = open("{o}.{c}.fasta".format(o=output_prefix, c=c), 'w')
                        f_d[c].write(">{0}\n{1}\n".format(newid, r.seq))


    if faa_file is not None:
        f_d = {}  # cluster --> file handle
        # writer pbid master file
        with open(output_prefix + '.faa', 'w') as f:
            for r in SeqIO.parse(open(faa_file), 'fasta'):
                if r.id in per_pbid:
                    newid = "{pbid}|{gene}|{transcript}|{clusters}".format(\
                            pbid=r.id,
                            gene=per_pbid[r.id]['gene'],
                            transcript=per_pbid[r.id]['transcript'],
                            clusters=";".join(per_pbid[r.id]['clusters']))
                    f.write(">{0}\n{1}\n".format(newid, r.seq))
                    for c in per_pbid[r.id]['clusters']:
                        if c not in f_d:
                            f_d[c] = open("{o}.{c}.faa".format(o=output_prefix, c=c), 'w')
                        f_d[c].write(">{0}\n{1}\n".format(newid, r.seq))
        for handle in f_d.values(): handle.close()


    if gff_file is not None:
        f_d = {}  # cluster --> file handle
        # writer pbid master file
        with open(output_prefix + '.gff', 'w') as f:
            for r in collapseGFFReader(gff_file):
                if r.seqid in per_pbid:
                    newid = "{pbid}|{gene}|{transcript}|{clusters}".format(\
                            pbid=r.seqid,
                            gene=per_pbid[r.seqid]['gene'],
                            transcript=per_pbid[r.seqid]['transcript'],
                            clusters=";".join(per_pbid[r.seqid]['clusters']))
                    write_collapseGFF_format(f, r)
                    for c in per_pbid[r.seqid]['clusters']:
                        if c not in f_d:
                            f_d[c] = open("{o}.{c}.gff".format(o=output_prefix, c=c), 'w')
                        write_collapseGFF_format(f_d[c], r)
        for handle in f_d.values(): handle.close()

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser("De-duplicate FLNC reads per cluster")
    parser.add_argument("corrected_csv", help="Annotated, error-corrected FLNC CSV file")
    parser.add_argument("cluster_file", help="Short read barcode to cluster CSV file")
    parser.add_argument("--fasta", help="(Optional) Fasta file (IDs should be PB.X.Y)")
    parser.add_argument("--gff", help="(Optional) GFF file (IDs should be PB.X.Y)")
    parser.add_argument("--faa", help="(Optional) Faa file (IDs should be PB.X.Y)")

    args = parser.parse_args()

    if not os.path.exists(args.corrected_csv):
        print("Input file {0} does not exist! Abort!".format(args.corrected_csv), file=sys.stderr)
        sys.exit(-1)

    if not os.path.exists(args.cluster_file):
        print("Input file {0} does not exist! Abort!".format(args.cluster_file), file=sys.stderr)
        sys.exit(-1)

    cluster_info = {}
    reader = DictReader(open(args.cluster_file), delimiter='\t')
    if ('cell_barcode' not in reader.fieldnames) or ('cluster' not in reader.fieldnames):
        print("Cluster file {0} must contain 'cell_barcode' and 'cluster' fields!".format(args.cluster_file), file=sys.stderr)
        sys.exit(-1)
    for r in reader:
        if r['cluster']!='NA':
            cluster_info[r['cell_barcode_rev']] = r['cluster']

    output_prefix = args.corrected_csv[:args.corrected_csv.rfind('.')] + '.dedup'

    main(args.corrected_csv, cluster_info, output_prefix, args.fasta, args.gff, args.faa)
