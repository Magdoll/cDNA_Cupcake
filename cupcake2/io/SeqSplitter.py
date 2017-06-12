#!/usr/bin/env python
"""Define Class `FastaSplitter` which splits a fasta file into
smaller files each containing `reads_per_split` reads."""

import os
import os.path as op
import sys

from Bio import SeqIO

from pbtranscript.Utils import mkdir



class FaFqSplitter(object):
    """
    Splits a fasta/fastq file into
    smaller chunks with a given prefix.
    """

    def __init__(self, input_fa_or_fq, reads_per_split, out_dir, out_format, is_fq):
        self.input_fa_or_fq = input_fa_or_fq
        self.is_fq = is_fq
        self.out_dir = out_dir
        self.reads_per_split = reads_per_split  # Number of reads per split
        self.out_format = out_format
        self.out_fns = None
        mkdir(self.out_dir)

    def __str__(self):
        if self.out_fns is None or len(self.out_fns) == 0:
            return "{input_fasta} ".format(input_fasta=self.input_fa_or_fq) + \
                "will be splitted into files each has " + \
                "{n} reads.".format(n=self.reads_per_split)
        else:
            return "{input_fasta} has been splitted into ".\
                   format(input_fasta=self.input_fa_or_fq) + \
                   "{m} files each has {n} reads:\n".\
                   format(m=len(self.out_fns),
                          n=self.reads_per_split) + ";".join(self.out_fns)

    def _out_fn(self, split_index):
        """Return name of the `split_index`-th splitted file."""
        if split_index > 999:
            raise ValueError("Too many splitted files to generate: number " +
                             "of splitted files exceed 1000.")
        name = self.out_format.format(split_index)
        return op.join(self.out_dir, name)

    def split(self, reads_in_first_split=None):
        """Split `input_fasta` into smaller files each containing
        `reads_per_split` reads. Return splitted fasta."""
        split_index = 0
        self.out_fns = []
        writer = open(self._out_fn(split_index), 'w')
        self.out_fns.append(self._out_fn(split_index))
        if reads_in_first_split is None:
            reads_in_first_split = self.reads_per_split

        io_format = 'fastq' if self.is_fq else 'fasta'
        reader = SeqIO.parse(open(self.input_fa_or_fq), io_format)
        for ridx, r in enumerate(reader):
            if ((split_index == 0 and ridx == reads_in_first_split) or
                    (split_index > 0 and ridx % self.reads_per_split == 0)) \
                and ridx != 0:
                split_index += 1
                writer.close()
                writer = open(self._out_fn(split_index), 'w')
                self.out_fns.append(self._out_fn(split_index))
            SeqIO.write(r, writer, io_format)
        writer.close()
        return list(self.out_fns)

    def rmOutFNs(self):
        """Remove splitted files."""
        for f in self.out_fns:
            os.remove(f)
        self.out_fns = []


def splitFaFq(input_fa_or_fq, reads_per_split, out_dir, out_format, is_fq, reads_in_first_split=None):
    """
    Split input_fasta into small fasta files each containing at most
    reads_per_split reads. All splitted fasta files will be placed under
    out_dir with out_prefix. Return paths to splitted files in a list.
    """
    obj = FaFqSplitter(input_fa_or_fq=input_fa_or_fq,
                        reads_per_split=reads_per_split,
                        out_dir=out_dir, out_format=out_format,
                        is_fq=is_fq)
    return obj.split(reads_in_first_split=reads_in_first_split)


def get_args():
    """Get arguments."""
    import argparse
    parser = argparse.ArgumentParser(
        description="Split a fasta/fastq file into smaller chunks.")
    parser.add_argument("input_fa_or_fq",
                        type=str,
                        help="Input fasta/fastq to be splitted.")
    parser.add_argument("reads_per_split",
                        type=int,
                        help="Reads per split.")
    parser.add_argument("out_dir",
                        type=str,
                        help="Output directory.")
    parser.add_argument("out_format",
                        type=str,
                        help="Output files format. (ex: \"input_split.{0:03d}.fasta\"")
    parser.add_argument("--is_fq", default=False, action="store_true", help="Is fastq (default: False)")
    this_args = parser.parse_args()
    return this_args


def main():
    """Main function, split a fasta into smaller chunks."""
    import logging
    from pbtranscript.__init__ import get_version
    log = logging.getLogger(__name__)
    args = get_args()
    from pbtranscript.Utils import setup_log
    setup_log(alog=log, level=logging.DEBUG)
    log.info("Running {f} v{v}.".format(f=op.basename(__file__),
                                        v=get_version()))

    splitFaFq(input_fa_or_fq=args.input_fa_or_fq,
               reads_per_split=args.reads_per_split,
               out_dir=args.out_dir,
               out_format=args.out_format,
               is_fq=args.is_fq)

if __name__ == "__main__":
    sys.exit(main())