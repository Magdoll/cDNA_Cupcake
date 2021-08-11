__author__ = "etseng@pacb.com"
"""
Tagging BAM files with phasing info 
"""

import pysam
from csv import DictReader


def main(read_bam, hap_info, output_bam, celltype_file=None):
    d = {}
    celltype_info = {}
    #for r in DictReader(open('phased.partial.cleaned.hap_info.txt'),delimiter=','):
    for r in DictReader(open(hap_info), delimiter=','):
        d[r['id']] = r['hap_postclean']

    if celltype_file is not None:
        for r in DictReader(open(celltype_file), delimiter=','):
            if r['id'] in d: celltype_info[r['id']] = r['Celltype'].replace(' ','_').replace(':','_')

    reader = pysam.AlignmentFile(read_bam, 'rb', check_sq=False)
    f2 = pysam.AlignmentFile(output_bam, 'wb', header=reader.header)
    for r in reader:
        d2 = r.to_dict()
        if r.qname in d:
            d2['tags'].append('RG:Z:' + str(d[r.qname]))
            if r.qname in celltype_info: d2['tags'].append('XC:Z:' + str(celltype_info[r.qname]))
        else: d2['tags'].append('RG:Z:NA')
        x = pysam.AlignedSegment.from_dict(d2, r.header)
        f2.write(x)

    f2.close()

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser("Tagging BAM files with phasing info")
    parser.add_argument("read_bam", help="Aligned BAM file that be tagged")
    parser.add_argument("hap_info", help="Comma-delimited hap info CSV, must have column 'id' and 'hap_postclean'")
    parser.add_argument("output_bam", help="Output tagged BAM filename")
    parser.add_argument("--celltype", default=None, help="[Optional] Comma-delimited celltype info CSV, must have column 'id' and 'Celltype'")

    args = parser.parse_args()
    main(args.read_bam, args.hap_info, args.output_bam, args.celltype)
