import os, re, sys
import subprocess

try:
    import vcf
except ImportError:
    print("Cannot import vcf! Please install pyvcf!", file=sys.stderr)
    sys.exit(-1)

from phasing.io import SAMMPileUpReader as sam
from phasing.io import MPileUpVariantCaller as VC
from phasing.io import VariantPhaser


MIN_COVERAGE = 10     # minimum number of FL reads for a gene to do SNP calling and phasing
ERR_SUB = 0.005
PVAL_CUTOFF = 0.1

def parse_user_input():
    from argparse import ArgumentParser
    parser = ArgumentParser(
            description = "A pipeline for aligning sequence data on a slurm cluster"
            )
    parser.add_argument('-a', '--assembly',
                        help="The mag assembly file in fasta format",
                        required=True, type=str
                        )
    parser.add_argument('-b', '--bamfile',
                        help="Aligned reads in bam file format [full path needed!]",
                        required=True, type=str
                        )
    parser.add_argument('-o', '--output',
                        help="output prefix",
                        required=True, type=str
                        )
    parser.add_argument('-g', '--genes',
                        help='SCG gene bed file',
                        required =True, type=str
                        )
    parser.add_argument('-p', '--pval_cutoff',
                        help="P value cutoff for variant calls",
                        default=PVAL_CUTOFF, type=float
                        )
    parser.add_argument("--bhFDR", default=None,
                        type=float,
                        help="FDR to be used for the Benjamini–Hochberg correction. Default: None (not used).")


    return parser.parse_args(), parser

def main(args, parser):
    args = parser.parse_args()

    if args.bhFDR is not None:
        print("--bhFDR {0} is given! Will be using Benjamini–Hochberg correction insteaad. --pval_cutoff is ignored.".format(args.bhFDR))


    # remove potential past run output
    past_files = [args.output+'.NO_SNPS_FOUND',
             args.output+'.NO_HAPS_FOUND',
             args.output+'.snps',
             args.output+'.log',
             args.output+'.human_readable.txt',
             args.output+'.vcf',
             args.output+'.cleaned.human_readable.txt',
             args.output+'.cleaned.vcf']

    for file in past_files:
        if os.path.exists(file):
            os.remove(file)

    snpsfound = False
    # (0) generate pileups
    f_human1 = open(args.output + '.human_readable_by_pos.txt', 'w')
    f_human1.write("haplotype\thapIdx\tcontig\tpos\tvarIdx\tbase\tcount\n")
    f_human2 = open(args.output + '.human_readable_by_hap.txt', 'w')
    f_human2.write("haplotype\thapIdx\tcontig\tcount\n")
    f_human3 = open(args.output + '.human_readable_by_read.txt', 'w')
    f_human3.write("read_id\thaplotype\thapIdx\n")

    for mpileupFile, contig, start, end in elitePileups(args.bamfile, args.genes, args.assembly, args.output):
        # (1) read the mpileup and vall variants
        reader = sam.MPileUpReader(mpileupFile)
        recs = [r for r in reader]
        vc = VC.MagMPileUPVariant(recs, min_cov=MIN_COVERAGE, err_sub=ERR_SUB, expected_strand='+-',
                                  pval_cutoff=args.pval_cutoff,
                                  bhFDR=args.bhFDR)
        vc.call_variant()
        print(vc.variant)

        if len(vc.variant) != 0:
            snpsfound = True
        else:
            continue

        # we write SNPs with the bases separated by "/" not "|" becuz we haven't phased them yet
        with open(args.output + '.snps', 'a+') as f_snp:
            for pos, v in vc.variant.items():
                f_snp.write("{contig}\t{pos}\t{bases}\t{counts}\n".format(\
                    contig=contig,\
                    pos=pos+1,\
                    bases="/".join([b for (b,c) in v]),\
                    counts="/".join([str(c) for (b,c) in v])))

        # (2) for each CCS read, assign a haplotype (or discard if outlier)
        pp = VariantPhaser.MagVariantPhaser(vc)
        pp.phase_variant(args.bamfile, [contig, start, end], args.output, partial_ok=True)
        print(pp.haplotypes)
        pp.haplotypes.get_haplotype_vcf_assignment()
        pp.haplotypes.write_haplotype_to_humanreadable(contig, f_human1, f_human2, f_human3, pp.seq_hap_info)
        os.remove(mpileupFile)
    f_human1.close()
    f_human2.close()
    f_human3.close()

    if not snpsfound:
        os.system("touch {out}.NO_SNPS_FOUND".format(out=args.output))
        os.remove(args.output + '.human_readable.txt')
        print("No SNPs found. END.", file=sys.stderr)


def elitePileups(aligned_bam : str, gene_bed : str, assembly : str, outprefix : str) -> str:
    """

    :param aligned_bam:
    :param gene_bed: gene bed to extract for making pileup
    :param assembly:
    :param outprefix:
    :return:
    """
    for line in open(gene_bed):
        # contig_4047     8476    8850    contig_4047_5
        chrom, s, e, name = line.strip().split()

        outfile = "{p}.{c}_{s}_{e}.pileup".format(p=outprefix, c=chrom, s=s, e=e)
        cmd = "samtools mpileup -r {c}:{s}-{e} -f {asm} -s {bam} > {o}".format(\
            c=chrom, s=s, e=e, asm=assembly, bam=aligned_bam, o=outfile)
        if subprocess.check_call(cmd, shell=True)!=0:
            print("FAILED TO RUN CMD: {0}. Abort!".format(cmd))
            sys.exit(-1)
        yield outfile, chrom, int(s), int(e)


if __name__ == "__main__":
    args, parser = parse_user_input()
    main(args, parser)
