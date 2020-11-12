import os, sys
import vcf
import glob
from collections import defaultdict, Counter

def main(dirs, vcf_filename='phased.partial.vcf', output='IsoSeq_IsoPhase.vcf'):
    no_snp_found_filename = vcf_filename[:vcf_filename.rfind('.')] + '.NO_SNPS_FOUND'
    snps_by_chrom = defaultdict(lambda: [])

    samp_ft = vcf.model.make_calldata_tuple(['GT', 'HQ'])

    for d in dirs:
        filename = os.path.join(d, vcf_filename)
        if not os.path.exists(filename):
            if not os.path.exists(no_snp_found_filename):
                print("VCF file {0} does not exist. Skipping.".format(filename), file=sys.stderr)
            continue
        reader = vcf.VCFReader(open(filename))

        for r in reader:
            c = Counter() # genotype -> count
            for x in r.samples:
                if x.data.GT.count('|') == 0:
                    c[x.data.GT] += x.data.HQ
                else:
                    for i,gt in enumerate(x.data.GT.split('|')):
                        c[gt] += x.data.HQ[i]
            c_keys = c.keys()
            genotype = "|".join(str(k) for k in c_keys)
            counts = ",".join(str(c[k]) for k in c_keys)
            r.samples = []
            r.samples = [vcf.model._Call(r, 'SAMPLE', samp_ft(*[genotype, counts]))]
            snps_by_chrom[r.CHROM].append((r.POS, r))

    keys = list(snps_by_chrom.keys())
    keys.sort()

    reader.samples = ['SAMPLE']
    f = vcf.Writer(open(output, 'w'), reader)
    for k in keys:
        v = snps_by_chrom[k]
        v.sort(key=lambda x: x[0])
        for pos,rec in v:
            f.write_record(rec)

    f.close()
    print("Output written to:", output)

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("-d", "--dir", default='by_loci/', help="Directory containing the subdirs of IsoPhase (default: by_loci/)")
    parser.add_argument("-o", "--output", default='IsoSeq_IsoPhase.vcf', help="Output VCF filename (default: IsoSeq_IsoPhase.vcf)")
    parser.add_argument("--vcf", default='phased.partial.cleaned.vcf',
                        choices=['phased.partial.vcf', 'phased.partial.cleaned.vcf',
                                 'phased.nopartial.vcf', 'phased.nopartial.cleaned.vcf'],
                        help="VCF to use per directory (default: phased.partial.cleaned.vcf)")
    parser.add_argument("-s", "--select_dirs", default=None, help="Comma separate list of directories to tally - if this is used, <dir> is ignored")


    args = parser.parse_args()

    if args.select_dirs is None:
        dirs = glob.glob(args.dir + '/*size*')
    else:
        dirs = args.select_dirs.split(',')
        for d in dirs:
            if not os.path.isdir(d) or not os.path.exists(d):
                print("{0} is not a directory or does not exist! Abort!".format(d), file=sys.stderr)
                sys.exit(-1)

    main(dirs, args.vcf, args.output)