#!/usr/bin/env python
__author__ = 'etseng@pacb.com'

"""
Convert a SAM file into GFF3 format (https://uswest.ensembl.org/info/website/upload/gff3.html)
 that also includes additional information that GMAP GFF3 for the 'mRNA' type:
   --- coverage
   --- identity
   --- matches
   --- mismatches
   --- indels


ex: from GMAP

6   cow_hereford    mRNA    109018744   109018861   .   +   .   \
     ID=myid.mrna1;Name=myid;Parent=myid.path1;coverage=29.5;identity=98.3;matches=116;mismatches=2;indels=0;unknowns=0

"""

#!/usr/bin/env python
import os, sys
import subprocess
from collections import Counter
from math import floor
from BCBio import GFF as BCBio_GFF
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from cupcake.io.BioReaders import  GMAPSAMReader

def convert_sam_rec_to_gff3_rec(r, qid_index_dict=None):
    """
    :param r: GMAPSAMRecord record
	:param qid_seen: list of qIDs processed so far -- if redundant, we have to put a unique suffix
    :return SeqRecord ready to be written as GFF3
    """
    if r.sID == '*':
        print >> sys.stderr, "Skipping {0} because unmapped.".format(r.qID)
        return None
    seq = Seq("AAAA")  # DO NOT CARE since sequence is not written in GFF3
    rec = SeqRecord(seq, r.sID)
    strand = 1 if r.flag.strand == '+' else -1

    indels = r.num_ins+r.num_del
    mismatches = r.num_nonmatches
    matches = r.num_mat_or_sub - r.num_nonmatches

    if qid_index_dict is not None:
        if r.qID in qid_index_dict: 
            qid_index_dict[r.qID] += 1
            r.qID += '_dup' + str(qid_index_dict[r.qID])
        else: qid_index_dict[r.qID] += 1

    gene_qualifiers = {"source": "hg38", "ID": r.qID+'.gene', "Name": r.qID} # for gene record
    mRNA_qualifiers = {"source": "hg38", "ID": r.qID, "Name": r.qID, "Parent": r.qID+'.gene',
                       "coverage": "{0:.2f}".format(r.qCoverage*10**2) if r.qCoverage is not None else "NA",
                       "identity": "{0:.2f}".format(r.identity*10**2),
                       "matches": matches, "mismatches": mismatches, "indels": indels}

    # gene line, one per record
    top_feature = SeqFeature(FeatureLocation(r.sStart, r.sEnd), type="gene", strand=strand, qualifiers=gene_qualifiers)
    # mRNA line, one per record
    top_feature.sub_features = [SeqFeature(FeatureLocation(r.sStart, r.sEnd), type="mRNA", strand=strand, qualifiers=mRNA_qualifiers)]
    # exon lines, as many exons per record
    for i,e in enumerate(r.segments):
        exon_qual = {"source": "hg38", "ID": "{0}.exon{1}".format(r.qID,i+1), "Name": r.qID, "Parent": r.qID}
        top_feature.sub_features.append(SeqFeature(FeatureLocation(e.start, e.end), type="exon", strand=strand, qualifiers=exon_qual))
    rec.features = [top_feature]
    return rec

def convert_sam_to_gff3(sam_filename, output_gff3, q_dict=None):
    qid_index_dict = Counter()
    with open(output_gff3, 'w') as f:
        recs = [convert_sam_rec_to_gff3_rec(r0, qid_index_dict) for r0 in GMAPSAMReader(sam_filename, True, query_len_dict=q_dict)]
        BCBio_GFF.write(filter(lambda x: x is not None, recs), f)

def main():
    from argparse import ArgumentParser

    parser = ArgumentParser("Convert SAM to GFF3 format using BCBio GFF")
    parser.add_argument("sam_filename")
    parser.add_argument("-i", "--input_fasta", default=None, help="(Optional) input fasta. If given, coverage will be calculated.")

    args = parser.parse_args()

    if not args.sam_filename.endswith('.sam'):
        print >> sys.stderr, "Only accepts files ending in .sam. Abort!"
        sys.exit(-1)

    prefix = args.sam_filename[:-4]
    output_gff3 = prefix + '.gff3'

    q_dict = None
    if args.input_fasta is not None:
        q_dict = dict((r.id, len(r.seq)) for r in SeqIO.parse(open(args.input_fasta), 'fasta'))

    with open(output_gff3, 'w') as f:
        recs = [convert_sam_rec_to_gff3_rec(r0) for r0 in GMAPSAMReader(args.sam_filename, True, query_len_dict=q_dict)]
        BCBio_GFF.write(filter(lambda x: x is not None, recs), f)


    print >> sys.stderr, "Output written to {0}.".format(output_gff3)

if __name__ == "__main__":
    main()
