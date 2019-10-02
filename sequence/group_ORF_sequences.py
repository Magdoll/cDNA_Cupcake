#!/usr/bin/env python
"""
Given a FAA file, group identical ORFs.

INPUT: a single FAA file
OUTPUT: a de-duped FAA file with a companion "group" file

The de-duped FAA format:

>ORFgroup_<PB.X>_<index>  if it is in PB.X.Y format

>ORFgroup_<index>  if not PB.X.Y format

The group file format:

ORFgroup_<index> \t comma-sep list of IDs with the same ORF
"""

import os, re, sys
from collections import defaultdict, OrderedDict, Counter
from Bio import SeqIO


rex_pbid = re.compile('(PB.\d+).(\d+)')
def dedup_ORFs(faa_filename, output_prefix, is_pbid):

    seq_dict = OrderedDict()  # ORF seq --> list of IDs

    pbid_counter = Counter()  # PB.X  --> counter  (not used if is_pbid is False)

    for r in SeqIO.parse(open(faa_filename), 'fasta'):
        s = str(r.seq).upper()
        if s not in seq_dict: seq_dict[s] = []
        seq_dict[s].append(r.id)

    f1 = open(output_prefix+'.faa', 'w')
    f2 = open(output_prefix+'.group.txt', 'w')

    for i,s in enumerate(seq_dict):
        newid = None
        if is_pbid:
            m = rex_pbid.match(seq_dict[s][0])  # we will just take the first member and use the PB.X.Y
            if m is None:
                print("WARNING: seqid {0} is not in PB.X.Y format!".format(seq_dict[s][0]), file=sys.stderr)
            else:
                pb_x = m.group(1)   # ex: PB.10
                pbid_counter[pb_x] += 1
                newid = "ORFgroup_{0}_{1}".format(pb_x, pbid_counter[pb_x])
        if newid is None:
            newid = "ORFgroup_{0}".format(i+1)
        f1.write(">{0}\n{1}\n".format(newid, s))
        f2.write("{0}\t{1}\n".format(newid, ",".join(seq_dict[s])))

    f1.close()
    f2.close()

    print("Output written to: {0},{1}".format(f1.name, f2.name), file=sys.stderr)


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser("De-duplicate ORF FAA file.")
    parser.add_argument("input_faa", help="Input FAA filename")
    parser.add_argument("output_prefix", help="Output prefix")
    parser.add_argument("--is_pbid", action="store_true", default=False, help="FAA IDs are in PB.X.Y format (default: off)")

    args = parser.parse_args()

    dedup_ORFs(args.input_faa, args.output_prefix, args.is_pbid)