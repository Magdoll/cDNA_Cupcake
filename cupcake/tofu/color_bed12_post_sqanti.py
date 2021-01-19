#!/home/UNIXHOME/etseng/anacondaPy37/envs/anaCogentPy37/bin/python
__author__ = 'etseng@pacb.com'

import pdb
import os, sys, math
from csv import DictReader
from collections import defaultdict, Counter

"""
Based on the script by Gloria Sheynkman for creating a BED12 file where isoforms are shaded by abundances

Read a SQANTI(3) classification output file that contains the FL abundance column (will also look at FL.<sample>)
Outputs a BED12 file
"""

# rgb scaling
RGB_SCALE = ['0,0,0', '26,0,0', '51,0,0', '77,0,0', '102,0,0',
             '128,0,0', '153,0,0', '179,0,0', '204,0,0', '230,0,0',
             '255,0,0', '255,26,26', '255,51,51', '255,77,77', '255,102,102',
             '255,128,128', '255,153,153', '255,179,179', '255,204,204', '255,230,230']
NUM_RGB = len(RGB_SCALE)

def shade_isoforms_for_gene_group(records, bed_info, bed_writers, ok_to_ignore=False):
    """

    :param isolist: list of SQANTI3 classification records which all come from same gene
    :return:
    """
    cpm_fields = bed_writers.keys()
    max_cpm_dict = {}
    for k in cpm_fields:
        max_cpm_dict[k] = max(r[k] for r in records)
    for r in records:
        if r['isoform'] not in bed_info and ok_to_ignore: continue
        info = bed_info[r['isoform']]
        for k in cpm_fields:
            if r[k] > 0:
                rgb_index = min(NUM_RGB-1, math.ceil(math.log2(max_cpm_dict[k]/r[k])*3))
            else:
                rgb_index = NUM_RGB-1
            # 0: chrom, 1: start, 2: end <-- no need to change
            # 3: name <-- change to gene|PB.X.Y|CPM
            info[3] = r['associated_gene'] + '|' + r['isoform'] + '|' + str(int(r[k]))
            # 4: score, 5: strand <-- no need to change
            #bed_writers[k].write(r['cdsStart'] + '\t')
            #bed_writers[k].write(r['cdsEnd'] + '\t')
            # 8: RGB color
            info[8] = RGB_SCALE[rgb_index]
            bed_writers[k].write("\t".join(info) + '\n')


def shaded_bed12_post_sqanti(sqanti_class_filename, input_bed12, output_prefix, FL_fieldnames=['FL'], ok_to_ignore=False):

    # read input BED12 file into dict
    bed_info = {}   # isoform --> bed record
    for line in open(input_bed12):
        raw = line.strip().split()
        bed_info[raw[3]] = raw

    CPM_fieldnames = {} # CPM -> FL field names
    for k in FL_fieldnames:
        assert k.startswith('FL')
        if k == 'FL':
            CPM_fieldnames["CPM"] = "FL"
        else:
            assert k.startswith('FL.')
            CPM_fieldnames["CPM."+k[3:]] = k

    # group SQANTI3 classification file by `associated_gene`
    records_by_gene = defaultdict(lambda: [])
    total_fl_count_dict = Counter()
    for r in DictReader(open(sqanti_class_filename), delimiter='\t'):
        records_by_gene[r['associated_gene']].append(r)
        for cpm_k, fl_k in CPM_fieldnames.items():
            total_fl_count_dict[cpm_k] += int(r[fl_k]) if r[fl_k]!='NA' else 0

    for cpm_k in total_fl_count_dict:
        if total_fl_count_dict[cpm_k] == 0:
            print("No counts observed in column `{0}`. Ignore!".format(CPM_fieldnames[cpm_k]), file=sys.stderr)
            del CPM_fieldnames[cpm_k]

    print("Generating count RGB for columns:".format(",".join(CPM_fieldnames.keys())))
    bed_writers = {}
    for cpm_k in CPM_fieldnames:
        outfile = output_prefix + '.' + cpm_k + '.bed12'
        print("Writing output to {0}....".format(outfile), file=sys.stdout)
        bed_writers[cpm_k] = open(outfile, 'w')
        bed_writers[cpm_k].write("track name=PacBioColored itemRgb=On\n")

    # calculate FL CPM
    for gene, records in records_by_gene.items():
        for r in records:
            for cpm_k, fl_k in CPM_fieldnames.items():
                r[cpm_k] = (int(r[fl_k]) if r[fl_k]!='NA' else 0) * (10**6) / total_fl_count_dict[cpm_k]
        shade_isoforms_for_gene_group(records, bed_info, bed_writers, ok_to_ignore)


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("class_filename", help='SQANTI(3) classification filename')
    parser.add_argument("bed12_filename", help='Input BED12 filename (converted from same SQANTI3 input GTF)')
    parser.add_argument("output_prefix")
    parser.add_argument("--ok_to_ignore", action="store_true", help="OK to ignore entries missing in bed file (default: off)")

    args = parser.parse_args()

    fl_fieldnames = []
    reader = DictReader(open(args.class_filename), delimiter='\t')
    for x in reader.fieldnames:
        if x.startswith('FL'):
            fl_fieldnames.append(x)

    if len(fl_fieldnames)==0:
        print("Expected column(s) 'FL' or 'FL.<sample>'! None found. Abort!")
        sys.exit(-1)

    shaded_bed12_post_sqanti(args.class_filename, args.bed12_filename, args.output_prefix, fl_fieldnames, args.ok_to_ignore)
