#!/usr/bin/env python
__author__ = 'etseng@pacb.com'

import os, sys
from cupcake.io.SeqReaders import LazyFastaReader
from phasing.io.MummerSNPReader import write_snp_to_vcf

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser("Process one or more .snps_files files from dnadiff to VCF format.")
    parser.add_argument("snps_filename", help="Filename containing list of .snps_files to process.")
    parser.add_argument("genome_filename", help="Genome fasta. Chromosome IDs must agree with .snps_files files!")

    args = parser.parse_args()

    snps_filename = args.snps_filename
    genome_filename = args.genome_filename

    snps_files = []
    # sanity checking of input files
    for line in open(snps_filename):
        filename = line.strip()
        if not filename.endswith('.snps'):
            print("Input files listed in {0} must end with .snps_files!".format(snps_filename), file=sys.stderr)
            sys.exit(-1)
        if not os.path.exists(filename):
            print("{0} does not exist! Abort.".format(filename), file=sys.stderr)
            sys.exit(-1)
        snps_files.append(filename)

    if not os.path.exists(genome_filename):
        print("Genome file {0} does not exist!".format(genome_filename), file=sys.stderr)

    print("Reading genome file {0}....".format(genome_filename), file=sys.stderr)
    genome_d = LazyFastaReader(genome_filename)

    # quick checking if the genome chromosomes have the |arrow|arrow style suffix, if they do, process it
    keys = list(genome_d.keys())
    for k in keys:
        k2 = k.split('|')[0]
        if k2!=k and k2 not in keys:
            genome_d.d[k2] = genome_d.d[k]
            print("Detected | string in chromosome ID, stripping {0} to {1}....".format(k, k2), file=sys.stderr)
    print("Finished reading genome.", file=sys.stderr)

    for snp_file in snps_files:
        assert snp_file.endswith('.snps')
        if os.stat(snp_file).st_size==0:
            print("Skipping {0} because empty file.".format(snp_file), file=sys.stderr)
            continue
        vcf_file = snp_file[:-5] + '.vcf'
        print("Processing {0} --> {1}".format(snp_file, vcf_file), file=sys.stderr)
        write_snp_to_vcf(snp_file, vcf_file, genome_filename, genome_d)


