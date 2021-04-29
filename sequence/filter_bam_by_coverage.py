#!/usr/bin/env python3
import pysam

def filter_bam_by_coverage(input_bam, output_bam, min_cov=0.9, filter_secondary=False, filter_supp=False):
    reader = pysam.AlignmentFile(input_bam, 'rb', check_sq=False)
    writer = pysam.AlignmentFile(output_bam, 'wb', template=reader)

    for r in reader:

        # calculate coverage, taking care of soft/hard clips
        # for cigar_stats():
        # The output order in the array is “MIDNSHP=X” followed by a field for the NM tag.
        # If the NM tag is not present, this field will always be 0.

        stats = r.get_cigar_stats()
        num_soft_hard_clipped = stats[0][4] + stats[0][5]
        cov = r.query_alignment_length * 1. / (r.query_alignment_length + num_soft_hard_clipped)
        if cov < min_cov:
            print("FILTER: {0} because coverage {1} too low.".format(r.qname, cov))
            continue
        if filter_secondary and r.is_secondary:
            print("FILTER: {0} because is secondary.".format(r.qname))
            continue
        if filter_supp and r.is_supplementary:
            print("FILTER: {0} because is supplementary.".format(r.qname))
            continue
        writer.write(r)
    reader.close()
    writer.close()
    print("Output written to:", output_bam)

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser("Filter BAM by coverage")
    parser.add_argument("input_bam")
    parser.add_argument("output_bam")
    parser.add_argument("-c", "--min_coverage", type=float, required=True, help="Minimum alignment coverage (between 0-1)")
    parser.add_argument("-s", "--filter_secondary", action="store_true", default=False, help="Filter secondary alignemnts (default: off)")
    parser.add_argument("-p", "--filter_supp", action="store_true", default=False, help="Filter supplementary alignemnts (default: off)")

    args = parser.parse_args()
    filter_bam_by_coverage(args.input_bam, args.output_bam, args.min_coverage, args.filter_secondary, args.filter_supp)
