#!/usr/bin/env python

__version__ = '1.0'

import os, sys
from Bio.Seq import Seq
from Bio import SeqIO
import coordinate_mapper as sp
import BioReaders


###### MODIFY FILENAME BELOW #######
#genome_file = 'hg38.fa'
#sam_file = 'alz.rep.fq.sam'
#sam_file = 'isoseq_flnc.fasta.sam'
#output_err_corrected_fasta = 'alz.rep.genome_corrected.fasta'
#output_err_corrected_fasta = 'isoseq_flnc.genome_corrected.fasta'
####################################

def err_correct(genome_file, sam_file, output_err_corrected_fasta, genome_dict=None):
    if genome_dict is None:
        genome_dict = {}
        for r in SeqIO.parse(open(genome_file), 'fasta'):
            genome_dict[r.name] = r
        print >> sys.stderr, "done reading", genome_file

    f = open(output_err_corrected_fasta, 'w')
    reader = BioReaders.GMAPSAMReader(sam_file, True)
    for r in reader:
        if r.sID == '*': continue
        seq = sp.consistute_genome_seq_from_exons(genome_dict, r.sID, r.segments, r.flag.strand)
        f.write(">{0}\n{1}\n".format(r.qID, seq))

    f.close()

    print >> sys.stderr, "output written to", output_err_corrected_fasta




if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser("Generate sequences using genome bases and SAM alignment file")
    parser.add_argument("genome_file",  help="Genome Fasta File")
    parser.add_argument("sam_file", help="GMAP SAM File")
    parser.add_argument("output_file", help="Output Fasta File")
    args = parser.parse_args()

    err_correct(args.genome_file, args.sam_file, args.output_file)