__author__ = 'etseng@pacb.com'

"""
For parsing the .snps_files results from running Mummer dnadiff into VCF
https://github.com/mummer4/mummer/blob/master/MANUAL.md#dnadiff

.snps_files file format:
(0) 1-based ref position
(1) base in ref
(2) base in query
(3) 1-based query position
(4) dist to nearest SNP
(5) ???
(6) length of ref
(7) length of query
(8) reading frame of ref (1 is forward, -1 reverse)
(9) reading frame of query
(10) name of ref
(11) name of query

IMPORTANT: Assumes ONLY snps_files and simple indels are in .snps_files. NO LARGE SVS!! This script does not handle them.

There are several diff between .snps_files format and VCF format (see Section 5 for examples)
https://samtools.github.io/hts-specs/VCFv4.2.pdf

for insertions, .snps_files shows ref pos but not the base. VCF expects the ref pos.

46060901        .       C       1679    0       1679    56369465        49753   1       1       000003F 000003F_064
46060901        .       C       1680    0       1680    56369465        49753   1       1       000003F 000003F_064

suppose ref@46060901 is "T"
this in VCF means: pos 46060901, ref: "T", alt: "TCC"

for deletions, .snps_files shows the ref pos but not the query base before it. VCF expects the ref pos before it.

46075044        A       .       15826   1       15826   56369465        49753   1       1       000003F 000003F_064
46075045        G       .       15826   1       15826   56369465        49753   1       1       000003F 000003F_064

suppose ref@46075043 is "T"
this in VCF means: pos 46075043, ref: "TAG", alt: "T"
"""

__VCF_EXAMPLE__ = \
"""
##fileformat=VCFv4.2
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT
20      1       .       G       A,T     .       PASS    AF=0.5       GT
"""

import os, sys
import vcf
from cupcake.io.SeqReaders import LazyFastaReader

class SNPRecord(object):
    def __init__(self, ref_pos, query_pos, ref_base, query_base, ref_len, query_len, ref_frame, query_frame, ref_name, query_name):
        """
        bases should be 0-based! when printed, will be 1-based.
        """
        self.ref_pos = ref_pos
        self.query_pos = query_pos
        self.ref_base = ref_base
        self.query_base = query_base
        self.ref_len = ref_len
        self.query_len = query_len
        self.ref_strand = '+' if int(ref_frame)==1 else '-'
        self.query_strand = '+' if int(query_frame)==1 else '-'
        self.ref_name = ref_name
        self.query_name = query_name

    def __str__(self):
        return """
        ref_pos: {0}:{1}
        query_pos: {2}:{3}
        ref_base: {4}
        query_base: {5}
        """.format(self.ref_name, self.ref_pos+1, self.query_name, self.query_pos+1, self.ref_base, self.query_base)

class SNPReader(object):
    def __init__(self, filename):
        self.filename = filename
        self.f = open(filename)

    def __iter__(self):
        return self

    def __next__(self):
        cur = self.f.tell()
        line = self.f.readline()
        if self.f.tell() == cur:
            raise StopIteration
        return self.parseLine(line)

    def parseLine(self, line):
        raw = line.strip().split('\t')
        if len(raw)!=12:
            raise Exception("Expected to have 12 cols in MUMMER SNP record \
            but saw only {0}, abort! Line was: {1}".format(len(raw), line))
        return SNPRecord(ref_pos=int(raw[0])-1,
                         query_pos=int(raw[3])-1,
                         ref_base = raw[1],
                         query_base = raw[2],
                         ref_len = int(raw[6]),
                         query_len = int(raw[7]),
                         ref_frame = int(raw[8]),
                         query_frame = int(raw[9]),
                         ref_name = raw[10],
                         query_name = raw[11])


def write_snp_to_vcf(snp_filename, vcf_filename, genome_filename, genome_d=None):
    # read the genome is genome_d is not given
    if genome_d is None:
        genome_d = LazyFastaReader(genome_filename)

    # read the first SNP record so we know the query name
    snp_reader = SNPReader(snp_filename)
    snp_rec = next(snp_reader)
    sample_name = snp_rec.query_name
    cur_recs = [snp_rec]
    genome_rec = genome_d[snp_rec.ref_name]

    with open('template.vcf', 'w') as f:
        f.write(__VCF_EXAMPLE__ + '\n')
    reader = vcf.VCFReader(open('template.vcf'))
    reader.samples = [sample_name]
    f_vcf = vcf.Writer(open(vcf_filename, 'w'), reader)

    for r1 in snp_reader:
        if r1.ref_pos == cur_recs[-1].ref_pos:  # multi-nt insertion, keep recording
            cur_recs.append(r1)
        elif r1.query_base == '.' and cur_recs[-1].query_base == '.': # multi-nt deletion, keep recording
            cur_recs.append(r1)
        else: # time to write out the current set of records
            # multiple records mean it could be:
            # 1. multi-nucleotide insertions
            # 2. multi-nucleotide deletions

            if len(cur_recs) == 1 and cur_recs[0].ref_base!='.' and cur_recs[0].query_base!='.': # just a SNP record
                pos = cur_recs[0].ref_pos
                ref_base = cur_recs[0].ref_base
                alt_base = cur_recs[0].query_base
            elif cur_recs[0].ref_base=='.':
                # is a single or multi-nt insertions, must retrieve ref base from genome
                # ex: in out.snps_files it is . --> ATG
                # in VCF it should be T --> TATG (meaning insertion of ATG)
                pos = cur_recs[0].ref_pos
                ref_base = genome_rec[cur_recs[0].ref_pos]
                alt_base = ref_base + "".join(r.query_base for r in cur_recs)
            else:
                # is a single multi-nt deletions, we need to get one more ref base before the first deletion
                # ex: in out.snps_files it is GGG --> deletion
                # in VCF it should be TGGG --> T (meaning deletion of GGG)
                pos = cur_recs[0].ref_pos-1
                ref_base_prev = genome_rec[pos]
                ref_base = ref_base_prev + "".join(r.ref_base for r in cur_recs)
                alt_base = ref_base_prev

            rec = vcf.model._Record(CHROM=snp_rec.ref_name,
                                POS=pos+1,
                                ID='.',
                                REF=ref_base,
                                ALT=[vcf.model._Substitution(alt_base)],
                                QUAL='.', FILTER='PASS',
                                INFO={'AF':0.5},
                                FORMAT="GT",
                                sample_indexes=None)
            samp_ft = vcf.model.make_calldata_tuple(['GT'])
            rec.samples.append(vcf.model._Call(rec, sample_name, samp_ft(*["0|1"])))
            f_vcf.write_record(rec)
            if r1.ref_name != cur_recs[0].ref_name:
                genome_rec = genome_d[r1.ref_name]
            cur_recs = [r1]
    f_vcf.close()


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
        vcf_file = snp_file[:-5] + '.vcf'
        print("Processing {0} --> {1}".format(snp_file, vcf_file), file=sys.stderr)
        write_snp_to_vcf(snp_file, vcf_file, genome_filename, genome_d)


