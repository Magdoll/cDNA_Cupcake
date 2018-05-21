#!/usr/bin/env python

import os, sys
from Bio import SeqIO

"""
Helper script for taking:
 -- ccs.fastq
 -- isoseq_flnc.fasta
 -- (optional) isoseq_nfl.fasta

 and converting it to: isoseq_flnc.fastq and isoseq_nfl.fasta
"""

# ex: m150803_002149_42161_c100745121910000001823165807071563_s1_p0/14/31_1114_CCS
def makerec(r, s, e):
    if s < e:
        return r[s:e]
    else:
        r2 = r[e:s].reverse_complement()
        r2.id = r.id
        r2.description = r.description
        return r2

def main():
    if not os.path.exists("isoseq_flnc.fasta"):
        print >> sys.stderr, "isoseq_flnc.fasta not found! Abort."
        sys.exit(-1)

    if not os.path.exists("ccs.fastq"):
        print >> sys.stderr, "ccs.fastq not found! Abort."
        sys.exit(-1)

    if not os.path.exists("isoseq_nfl.fasta"):
        print >> sys.stderr, "WARNING: isoseq_nfl.fasta is not found. Skipping."
        nfl = {}
    else:
        nfl = dict((r.id[:r.id.rfind('/')],r.description) for r in SeqIO.parse(open('isoseq_nfl.fasta'),'fasta'))

    flnc = dict((r.id[:r.id.rfind('/')],r.description) for r in SeqIO.parse(open('isoseq_flnc.fasta'),'fasta'))
    f1 = open('isoseq_flnc.fastq', 'w')
    if len(nfl) > 0:
        f2 = open('isoseq_nfl.fastq', 'w')


    for r in SeqIO.parse(open('ccs.fastq'),'fastq'):
        p = r.id[:r.id.rfind('/')]
        if p in flnc:
            r.description = flnc[p]
            r.id = r.description.split()[0]
            s, e, junk = r.id.split('/')[-1].split('_')
            SeqIO.write(makerec(r, int(s), int(e)), f1, 'fastq')
            del flnc[p]
        elif p in nfl:
            r.description = nfl[p]
            r.id = r.description.split()[0]
            s, e, junk = r.id.split('/')[-1].split('_')
            SeqIO.write(makerec(r, int(s), int(e)), f2, 'fastq')
            del nfl[p]

    f1.close()
    print >> sys.stderr, "Output written to: isoseq_flnc.fastq"
    if len(nfl) > 0:
        print >> sys.stderr, "Output written to: isoseq_nfl.fastq"
        f2.close()


if __name__ == "__main__":
    main()





