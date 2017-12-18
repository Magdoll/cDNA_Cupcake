

__author__ = 'etseng@pacb.com'

import os
import sys
import vcf
from Bio import SeqIO
import phasing.io.SAMMPileUpReader as sp
import phasing.io.MPileUpVariantCaller as VC
from phasing.io import VariantPhaser
from phasing.io import VariantPhaseCleaner

MIN_COVERAGE = 10
ERR_SUB = 0.005
MAX_DIFF_ALLOWED = 3  # maximum difference in bases allowed for two haplotype strings
PVAL_CUTOFF = 0.01

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("fasta_filename")
parser.add_argument("sam_filename")
parser.add_argument("mpileup_filename")
parser.add_argument("read_stat")
parser.add_argument("mapping_filename")
parser.add_argument("-o", "--output_prefix", required=True)
parser.add_argument("--strand", choices=['+', '-'], required=True)
parser.add_argument("--partial_ok", default=False, action="store_true")
parser.add_argument("-p", "--pval_cutoff", default=PVAL_CUTOFF, type=float)
parser.add_argument("-n", "--ploidity", type=int, default=2)
#parser.add_argument("-e", "--err_sub", default=ERR_SUB, type=float, help="Estimated substitution error rate (default: 0.005)")


args = parser.parse_args()


# (1) read the mpileup and vall variants
reader = sp.MPileUpReader(args.mpileup_filename)
recs = [r for r in reader]
vc = VC.MPileUPVariant(recs, min_cov=MIN_COVERAGE, err_sub=ERR_SUB, expected_strand=args.strand, pval_cutoff=args.pval_cutoff)
vc.call_variant()
print vc.variant

if len(vc.variant) == 0:
    os.system("touch NO_SNPS_FOUND")
    print >> sys.stderr, "No SNPs found. END."
    sys.exit(0)

# (2) for each CCS read, assign a haplotype (or discard if outlier)
pp = VariantPhaser.VariantPhaser(vc)
pp.phase_variant(args.sam_filename, args.fasta_filename, args.output_prefix, partial_ok=args.partial_ok)
pp.haplotypes
pp.haplotypes.get_haplotype_vcf_assignment()

# (3) phase isoforms
seqids = set([r.id for r in SeqIO.parse(open(args.fasta_filename), 'fasta')])
isoform_tally = VariantPhaser.phase_isoforms(args.read_stat, seqids, pp)
if len(isoform_tally) == 0:
    os.system("touch NO_HAPS_FOUND")
    print >> sys.stderr, "No good haps found. END."
    sys.exit(0)
pp.haplotypes.write_haplotype_to_vcf(args.mapping_filename, isoform_tally, args.output_prefix)

# (4) clean isoforms
hap_count = VariantPhaseCleaner.make_haplotype_counts(isoform_tally)

# --- old, obsolete ----
#G, partial_haps = VariantPhaseCleaner.make_haplotype_graph_nonpartial_only(pp.haplotypes.haplotypes, ERR_SUB, MAX_DIFF_ALLOWED)
#print "G nodes:", G.nodes(), "G edges:", G.edges()
#m, new_hap, new_isoform_tally = VariantPhaseCleaner.error_correct_haplotypes(G, partial_haps, hap_count, pp.haplotypes, isoform_tally)
diff_arr, hap_count_ordered = VariantPhaseCleaner.infer_haplotypes_via_min_diff(pp.haplotypes.haplotypes, hap_count, args.ploidity, MAX_DIFF_ALLOWED)
m, new_hap, new_isoform_tally = VariantPhaseCleaner.error_correct_haplotypes(pp.haplotypes, isoform_tally, diff_arr, hap_count_ordered)

new_hap.get_haplotype_vcf_assignment()
new_hap.write_haplotype_to_vcf(args.mapping_filename, new_isoform_tally, args.output_prefix+'.cleaned')