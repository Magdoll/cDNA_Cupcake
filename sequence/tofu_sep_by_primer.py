#!/usr/bin/env python

__version__ = '1.0'

# ##################################################################
#
# Separate input fasta/fastq sequence by primer (ToFU/Iso-Seq format)
#
# Each fasta/fastq sequence must have header description `primer=<pid>;`
#  to be properly parsed. <pid> is expected to an integer.
#
# example usage: tofu_sep_by_primer.py isoseq_flnc.fasta isoseq_flnc
#
# use `--just_count` option if you just want to get the primer counts.
#
# ##################################################################

import os, sys
from collections import defaultdict
from Bio import SeqIO

def type_fa_or_fq(file):
    file = file.upper()
    if file.endswith('.FA') or file.endswith('.FASTA'): return 'fasta'
    else: return 'fastq'

def sep_by_primer(file, output_prefix, just_count=False):
    """
    Separate an input fasta/fastq file by primer.
    Use the primer= information in the sequence header.
    ex: m151104_115017_42133_c100929862550000001823209905251602_s1_p0/18/1991_72_CCS strand=-;fiveseen=1;polyAseen=1;three
seen=1;fiveend=39;polyAend=1958;threeend=1989;primer=3;chimera=0

    Outputs will be written to <output_prefix>.primer_XXX.fasta/fastq
    """
    filetype = type_fa_or_fq(file)
    print "file type is:", filetype

    handle_dict = {} # pid --> file handle
    count_per_primer = defaultdict(lambda: 0) # pid -> count

    for r in SeqIO.parse(open(file), type_fa_or_fq(file)):
        h = None
        for stuff in r.description.split(';'):
            if stuff.startswith('primer='):
                pid = stuff.split('=')[1]
                count_per_primer[pid] += 1
                if not just_count:
                    if pid in handle_dict: h = handle_dict[pid]
                    else:
                        h = open("{0}.primer_{1}.{2}".format(output_prefix, pid, filetype), 'w')
                        handle_dict[pid] = h
                break
        if h is not None:
            SeqIO.write(r, h, filetype)
        elif not just_count:
            print >> sys.stderr, "Sequence {0} has no primer= metadata. Ignore!".format(r.description)

    if not just_count:
        for h in handle_dict.itervalues():
            print >> sys.stderr, "Primer file written:", h.name
            h.close()

    for pid, count in count_per_primer.iteritems():
        print "Count for primer {0}: {1}".format(pid, count)


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser("Separate fasta/fastq files by primer information")
    parser.add_argument("filename", help="Input fasta/fastq filename")
    parser.add_argument("output_prefix", help="Output file prefix")
    parser.add_argument("--just_count", default=False, action="store_true", help="Just provide primer counts, don't actually create the per-primer split files.")
    args = parser.parse_args()

    sep_by_primer(args.filename, args.output_prefix, args.just_count)


