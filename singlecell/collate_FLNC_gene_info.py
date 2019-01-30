
"""
Given:

   1. single cell UMI, BC file (ccs -> umi, bc)
   2. collapse group.txt  (ccs -> pbid)
   3. SQANTI classification (pbid -> transcript, isoform, category)
   4. optional ontarget file (pbid -> ontarget or not)

Output a collated infor file that is:

   <ccs>, <pbid>, <transcript>, <gene>, <category>, <ontarget Y|N|NA>, <UMI>, <BC>
"""

import os, sys
from csv import DictReader, DictWriter

def read_group_info(group_filename):
    """
    :return: dict of ccs -> pbid
    """
    d = {}
    for line in open(group_filename):
        pbid, members = line.strip().split('\t')
        for m in members.split(','):
            d[m] = pbid
    return d

def collate_gene_info(group_filename, csv_filename, class_filename, output_filename, ontarget_filename=None):
    """
    <id>, <pbid>, <transcript>, <gene>, <category>, <ontarget Y|N|NA>, <UMI>, <BC>
    """
    FIELDS = ['id', 'pbid', 'transcript', 'gene', 'category', 'ontarget', 'UMI', 'BC']

    group_info = read_group_info(group_filename)
    umi_bc_info = dict((r['id'], r) for r in DictReader(open(csv_filename), delimiter='\t'))
    sqanti_info = dict((r['isoform'], r) for r in DictReader(open(class_filename), delimiter='\t'))
    if ontarget_filename is not None:
        ontarget_info = dict((r['read_id'], r) for r in DictReader(open(ontarget_filename), delimiter='\t'))

    f = open(output_filename, 'w')
    writer = DictWriter(f, FIELDS, delimiter='\t')
    writer.writeheader()

    for ccs_id, pbid in group_info.iteritems():
        if pbid not in sqanti_info:
            print >> sys.stderr, "ignoring ID {0} cuz not in classification file.".format(pbid)
            continue
        rec = {'id': ccs_id, 'pbid': pbid}
        rec['category'] = sqanti_info[pbid]['structural_category']
        rec['transcript'] = sqanti_info[pbid]['associated_transcript']
        rec['gene'] = sqanti_info[pbid]['associated_gene']
        rec['UMI'] = umi_bc_info[ccs_id]['UMI']
        rec['BC'] = umi_bc_info[ccs_id]['BC']
        if ontarget_filename is None:
            rec['ontarget'] = 'NA'
        else:
            rec['ontarget'] = 'Y' if ontarget_info[pbid]['genes']!='' else 'N'

        writer.writerow(rec)

    f.close()


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("group_filename", help="Collapse .group.txt")
    parser.add_argument("csv_filename", help="Trimmed UMI/BC CSV info")
    parser.add_argument("class_filename", help="SQANTI classification.txt")
    parser.add_argument("output_filename", help="Output filename")
    parser.add_argument("-i", "--ontarget_filename", help="Optional on target information text")

    args = parser.parse_args()

    if os.path.exists(args.output_filename):
        print >> sys.stderr, "Output file {0} already exists. Abort!".format(args.output_filename)
        sys.exit(-1)

    if not os.path.exists(args.group_filename):
        print >> sys.stderr, "Group file {0} not found. Abort!".format(args.group_filename)
        sys.exit(-1)

    if not os.path.exists(args.csv_filename):
        print >> sys.stderr, "CSV file {0} not found. Abort!".format(args.csv_filename)
        sys.exit(-1)

    if not os.path.exists(args.class_filename):
        print >> sys.stderr, "Class file {0} not found. Abort!".format(args.class_filename)
        sys.exit(-1)

    collate_gene_info(args.group_filename, args.csv_filename, args.class_filename, args.output_filename, args.ontarget_filename)
