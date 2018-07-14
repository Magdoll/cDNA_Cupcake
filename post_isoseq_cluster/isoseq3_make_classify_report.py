__author__ = 'etseng@pacb.com'

"""
Helper script for making `classify_report.csv` from lima output for Iso-Seq.

lima must be run using parameters `--isoseq --dump-clips`.
lima clips look like:

>m54020_170625_150952/69664969/0_39 bq:85 bc:7
AAGCGTGGTATAACGCAGAGTACTCTGTATCTCTATGTG
>m54020_170625_150952/69664969/2692_2723 bq:80 bc:0
CCCCATGTAATCTGCGTTGATACCACGGCTT
>m54020_170625_150952/69664996/0_40 bq:80 bc:6
AAGACAGTGGTATCAACGCAGAGTACCATAATATCACTGT
>m54020_170625_150952/69664996/1058_1089 bq:100 bc:0
CCCCATGTACTCTGCGTTGATACCACTGCTT

classify_report.csv format:
id,strand,fivelen,threelen,polyAlen,insertlen,primer_index,primer
movie/zmw1,+,25,40,38,1000,0--1,Clontech_5p--Heart_3p
movie/zmw2,-',25,40,38,1200,0--2,Clontech_5p-Liver_3p
"""

import os, re, sys
from csv import DictWriter
from Bio import SeqIO
import pysam

clip_rex = re.compile('(\S+)/(\d+)_(\d+) bq:\d+ bc:(\d+)')

report_fields = ['id','strand','fivelen','threelen','polyAlen','insertlen','primer_index','primer']


def read_primer_fasta(primer_filename):
    primer_index_dict = {} # primer index --> (5p|3p, name)
    i = 0
    for r in SeqIO.parse(open(primer_filename), 'fasta'):
        if r.id.endswith('_5p'):
            primer_index_dict[i] = ('5p', r.id[:-3])
        elif r.id.endswith('_3p'):
            primer_index_dict[i] = ('3p', r.id[:-3])
        else:
            print >> sys.stderr, "Unrecognized primer ID format! Primes must end with _5p or _3p!"
            sys.exit(-1)
        i += 1
    return primer_index_dict


def make_classify_report_from_lima(clips_filename, primer_index_dict, flnc_bam=None):
    """
    clips format:

    """
    if flnc_bam is not None:
        flnc_len_dict = dict((r.qname, r.qlen) for r in pysam.Samfile(flnc_bam, check_sq=False))
    else:
        flnc_len_dict = None
        print >> sys.stderr, "WARNING: FLNC BAM not provided. `polyAlen` and `insertlen` fields will be `NA`."


    f = open('classify_report.csv', 'w')
    writer = DictWriter(f, fieldnames=report_fields, delimiter=',')
    writer.writeheader()
    first_of_pair_seen = False
    rec = {'id':None,'strand':None,'fivelen':None,'threelen':None,'polyAlen':None,'insertlen':None,
                   'primer_index':None,'primer':None}

    for r in SeqIO.parse(open(clips_filename),'fasta'):
        m = clip_rex.match(r.description)
        zmw = m.group(1) + '/ccs'
        s, e = int(m.group(2)), int(m.group(3))
        bc = int(m.group(4))

        if primer_index_dict[bc][0]=='5p':
            p5 = bc
            rec['fivelen'], start5, end5 = e-s, s, e
        else:
            assert primer_index_dict[bc][0]=='3p'
            p3 = bc
            rec['threelen'], start3, end3 = e-s, s, e

        if first_of_pair_seen:  # both pairs seen, write out and reset
            assert rec['id'] == zmw
            rec['strand'] = '+' if end5 < end3 else '-'
            if flnc_len_dict is None or zmw not in flnc_len_dict:
                rec['insertlen']='NA'
                rec['polyAlen']='NA'
            else:
                rec['insertlen'] = flnc_len_dict[zmw]
                if rec['strand'] == '+':
                    rec['polyAlen'] = start3 - end5 -rec['insertlen']
                else:
                    rec['polyAlen'] = start5 - end3 -rec['insertlen']

            rec['primer'] = "{0}--{1}".format(primer_index_dict[p5][1], primer_index_dict[p3][1])
            rec['primer_index'] = "{0}--{1}".format(p5, p3)
            writer.writerow(rec)
            #f.write("{id},{len5},{len3},{lenA},{lenI},{pn5}--{pn3},{p5}--{p3}\n".format(\
            #    id=zmw, len5=len5, len3=len3,
            #    lenA=lenA, lenI=lenI,
            #    pn5=primer_index_dict[p5][1], pn3=primer_index_dict[p3][1],
            #    p5=p5, p3=p3))

            # reset variables
            first_of_pair_seen = False
            rec = {'id':None,'strand':None, 'fivelen':None,'threelen':None,'polyAlen':None,'insertlen':None,
                   'primer_index':None,'primer':None}
            p5, p3, start5, end5, start3, end3 = None, None, None, None, None, None
        else: # first of the pair
            rec['id'] = zmw
            first_of_pair_seen = True
    f.close()
    print >> sys.stderr, "Classify report written to: {0}".format(f.name)


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser("Make classify_report.csv for Iso-Seq demultiplexed result using LIMA.")
    parser.add_argument("lima_clips", help="Clips output from lima")
    parser.add_argument("primer_fasta", help="Primer fasta file used to run LIMA")
    parser.add_argument("--flnc_bam", default=None, help="unpolished.flnc.bam (optional)")

    args = parser.parse_args()
    primer_index_dict = read_primer_fasta(args.primer_fasta)
    make_classify_report_from_lima(args.lima_clips, primer_index_dict, args.flnc_bam)