import os, sys, re
from csv import DictWriter
import pysam

report_fields = ['cluster_id', 'read_id', 'read_type']


def make_cluster_report_from_polished_bam(polished_bam):
    f = pysam.Samfile(polished_bam)
    h = open('cluster_report.csv', 'w')
    writer = DictWriter(h, fieldnames=report_fields, delimiter=',')
    writer.writeheader()
    for r in f:
        d = dict(r.tags)
        # transcript name format: transcript/0
        # zmw is of format <movie>/<zmw>/ccs, which is the ID for the FL read
        for zmw in d['im'].split(','):
            h.write("{0},{1},FL\n".format(r.qname, zmw))

    h.close()
    print >> sys.stderr, "Cluster report  written to: {0}".format(h.name)


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser("Make cluster_report.csv for Iso-Seq3 cluster output.")
    parser.add_argument("polished_bam", help="Polished bam file")

    args = parser.parse_args()

    if not os.path.exists(args.polished_bam):
        print >> sys.stderr, "{0} does not exist. Abort!".format(args.polished_bam)
        sys.exit(-1)

    make_cluster_report_from_polished_bam(args.polished_bam)