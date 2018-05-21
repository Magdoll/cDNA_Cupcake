#!/usr/bin/env python

"""
Settings should be provided via a config file containing

NAME=FMR1
SEQ_FILENAME=FMR1_L446_R503.fasta
BEFORE_SEQ_ID=FMR1_L446
AFTER_SEQ_ID=FMR1_R503
REPEAT=CGG;AGG;
SHOW_POS_FOR=AGG;
MIN_IDENTITY=0.9
MAX_MISSED_BEFORE_BP=10
MAX_MISSED_AFTER_BP=10


1. run BLASR to match against L/R (before/after) sequence
2. trim off the L/R before/after
3. identify location and number of repeats for each repeat(motif)

output:
(a) a summary file of <seqid>, <count_for_motif1>, <count_for_motif2>, etc...
(b) a fasta file where the repeats are masked out and only non-matches are shown

"""
import os, re, sys
from collections import namedtuple
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from cupcake.io.BLASRRecord import BLASRM4Reader
from cupcake.io.SeqReaders import LazyFastaReader
import repeat_counter as sp2

RepeatConfig = namedtuple('RepeatConfig', ['name', 'file', 'before_id', 'after_id',
                                           'motifs', 'min_identity', 'max_missed_before', 'max_missed_after',
                                           'show_pos_for'])

def read_config(config_filename):
    name, file, before_id, after_id, repeats = None, None, None, None, None
    min_iden = 0.9
    max_missed_before = 10
    max_missed_after = 10
    for line in open(config_filename):
        line = line.strip()
        if line.startswith('NAME='):
            name = line.split('=')[1]
        elif line.startswith('SEQ_FILENAME='):
            file = line.split('=')[1]
        elif line.startswith('BEFORE_SEQ_ID='):
            before_id = line.split('=')[1]
        elif line.startswith('AFTER_SEQ_ID='):
            after_id = line.split('=')[1]
        elif line.startswith('REPEAT='):
            repeats = line.split('=')[1].split(';')
            if len(repeats[-1]) == 0: repeats = repeats[:-1] # remove the last ;
        elif line.startswith('MIN_IDENTITY='):
            min_iden = float(line.split('=')[1])
            assert 0 < min_iden <= 1
        elif line.startswith('MAX_MISSED_BEFORE_BP='):
            max_missed_before = int(line.split('=')[1])
        elif line.startswith('MAX_MISSED_AFTER_BP='):
            max_missed_after = int(line.split('=')[1])
        elif line.startswith('SHOW_POS_FOR='):
            show_pos_for = line.split('=')[1].split(';')
            if len(show_pos_for[-1]) == 0: show_pos_for = show_pos_for[:-1] # remove last ;
        else:
            raise Exception, "Config line not recognized! {0} is not a valid config".format(line)

    if name is None:
        print >> sys.stderr, "NAME= is missing in config file. Abort!"
        sys.exit(-1)
    if file is None:
        print >> sys.stderr, "SEQ_FILENAME= is missing in config file. Abort!"
        sys.exit(-1)
    if before_id is None:
        print >> sys.stderr, "BEFORE_SEQ_ID= is missing in config file. Abort!"
        sys.exit(-1)
    if after_id is None:
        print >> sys.stderr, "AFTER_SEQ_ID= is missing in config file. Abort!"
        sys.exit(-1)
    if repeats is None:
        print >> sys.stderr, "REPEAT= is missing in config file. Abort!"
        sys.exit(-1)

    return RepeatConfig(name, file, before_id, after_id, repeats, min_iden, max_missed_before, max_missed_after, show_pos_for)


def iter_blasr_match(blasr_filename, before_id, after_id, min_identity=.90, max_missed_before_len=10, max_missed_after_len=10):
    """
    :param blasr_filename: must be a BLASR m4 header
    :return: (r1, r2) where both are same sequence, r1 is the one matching the 5' end, r2 matches 3' end
             note that it's possible it's matching on the - strand for the query sequence
    """
    reader = BLASRM4Reader(blasr_filename)
    it = iter(reader)
    prev = it.next()
    for r in it:
        if prev is not None and r.qID == prev.qID:
            r1 = prev
            r2 = r
            if r1.strand != r2.strand:
                print >> sys.stderr, "{0} pair matches have opp strands. Ignore.".format(r1.qID)
                continue
            elif r1.sID == r2.sID:
                print >> sys.stderr, "{0} has internal matches. Ignore".format(r1.qID)
                continue
            else:
                # confirm that one matches the before seq, one matches the after seq
                if r1.sID == before_id:
                    if r2.sID != after_id:
                        print >> sys.stderr, "{0} pair matches is not expected before/after. Ignore.".format(r1.qID)
                        continue
                elif r1.sID == after_id:
                    if r2.sID != before_id:
                        print >> sys.stderr, "{0} pair matches is not expected before/after. Ignore.".format(r1.qID)
                        continue
                    r1, r2 = r2, r1 # switch up so r1 is always the "before" and r2 the "after"
                else:
                    print >> sys.stderr, "{0} pair matches is not expected before/after. Ignore.".format(r1.qID)
                    continue
                # at this point, r1 is "before", r2 is "after" and strands do match
                # check identity
                if r1.identity/100. < min_identity or r2.identity/100. < min_identity:
                    print >> sys.stderr, "{0} pair matches identity too low. Ignore.".format(r1.qID)
                    continue
                # check that the max missed length is  under threshold
                missed_before = r1.sStart + (r1.sLength-r1.sEnd)
                missed_after = r2.sStart + (r2.sLength-r2.sEnd)
                if missed_before > max_missed_before_len or missed_after > max_missed_after_len:
                    print >> sys.stderr, "{0} pair matches coverage too low. Ignore.".format(r1.qID)
                    continue
                else: # finally! a good match with good identity and coverage
                    prev = None
                    yield (r1, r2)
        else:
            if prev is not None:
                # no paired match, ignore prev
                print >> sys.stderr, "{0} only has a single match. Ignore.".format(prev.qID)
            prev = r


def trim_seq(seqrec, r1, r2):
    """
    seqrec --- SeqRecord of the sequence to be trimmed
    r1 -- the "before" sequence match
    r2 -- the "after" sequence match

    return: a new SeqRecord with just the repeat region, before/after trimmed, and in correct orientation
    """
    assert r1.qID == r2.qID and r1.strand == r2.strand
    if r1.strand == '+':
        newrec = SeqRecord(Seq(seqrec[r1.qEnd:r2.qStart].seq))
    else:
        newrec = SeqRecord(Seq(seqrec[r2.qEnd:r1.qStart].seq)).reverse_complement()
    return newrec


def process_blasr_for_repeat(fasta_filename, blasr_filename, config):
    d = LazyFastaReader(fasta_filename)

    for (r1,r2) in iter_blasr_match(blasr_filename, config.before_id, config.after_id, \
                                    min_identity=config.min_identity,\
                                    max_missed_before_len=config.max_missed_before,\
                                    max_missed_after_len=config.max_missed_after):
        newrec = trim_seq(d[r1.qID], r1, r2)
        newseq = newrec.seq.tostring()
        intervals, tally = sp2.count_repeats(newseq, config.motifs)
        sp2.sanity_check(newseq, config.motifs, tally)
        masked_seq = sp2.mask_repeats_in_sequence(newseq, intervals)
        yield r1.qID, intervals, tally, masked_seq, newseq


def main(fasta_filename, blasr_filename, config_filename, output_prefix):
    config = read_config(config_filename)

    f1 = open(output_prefix + '.summary.txt', 'w')
    f1.write("seqid")
    for motif in config.motifs: f1.write("\tcount_{0}".format(motif))
    for motif in config.show_pos_for: f1.write("\tpos_{0}".format(motif))
    f1.write("\n")

    f2 = open(output_prefix + '.masked.fasta', 'w')
    f3 = open(output_prefix + '.trimmed.fasta', 'w')

    for seqid, intervals, tally, masked_seq, trimmed_seq in process_blasr_for_repeat(fasta_filename, blasr_filename, config):
        f1.write(seqid)
        for motif in config.motifs:
            f1.write("\t{0}".format(len(tally[motif])))
        for motif in config.show_pos_for:
            f1.write("\t{0}".format(",".join(map(lambda x: str(x+1), tally[motif]))))
        f1.write("\n")
        f2.write(">{0}\n{1}\n".format(seqid, masked_seq))
        f3.write(">{0}\n{1}\n".format(seqid, trimmed_seq))

    f1.close()
    f2.close()
    f3.close()

if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser("Identifying repeats in CCS data")
    parser.add_argument("fasta_filename")
    parser.add_argument("blasr_filename")
    parser.add_argument("config_filename")
    parser.add_argument("output_prefix")

    args = parser.parse_args()
    main(args.fasta_filename, args.blasr_filename, args.config_filename, args.output_prefix)