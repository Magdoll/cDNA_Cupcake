"""
Make FL (optionally, also with nFL) SAM alignment against the final isoform.

Uses the .read_stat.txt file to figure out which FL/nFL reads were assigned to which isoforms.
Then aligns the FL/nFL reads to the isoform sequence (ex: PB.1.1) using BLASR.

This is useful for manually examining the ICE algorithm is doing a good job and for debugging.

NOTE: in current ToFU/Iso-Seq, only the first 100 reads are used for Polishing since beyond that
      Quiver/Arrow don't get much improvement. However for this script we will align all the reads.

Requirement: BioPython, SeqReaders (through Cupcake), pbtranscript.ice (through pitchfork or installed SA)
"""

import os, sys
from csv import DictReader
from collections import defaultdict

from Bio import SeqIO
try:
    import SeqReaders        # get this through sequence/ dir added to PYTHONPATH
except ImportError:
    raise ImportError, "Unable to import SeqReaders! Please make sure cDNA_Cupcake/sequence is in your PYTHONPATH."

try:
    import pbtranscript.ice.IceUtils as sp   # must have SA3.2+ / pitchfork installed
    sp_version = 'NEW'
except ImportError: # try fall back to SA2.x / ToFU1
    try:
        import pbtools.pbtranscript.ice.IceUtils as sp
        sp_version = 'OLD'
    except ImportError:
        raise ImportError, "Unable to import pbtranscript.ice.IceUtils or pbtools.pbtranscript.ice.IceUtils! " \
                           "Please make sure you have SMRTAnalysis 2.x or SMRTAnalysis 3.x installed!"


def type_fa_or_fq(file):
    file = file.upper()
    if file.endswith('.FA') or file.endswith('.FASTA'): return 'fasta'
    elif file.endswith('.FQ') or file.endswith('.FASTQ'): return 'fastq'
    else:
        raise Exception, "File needs to end with .fasta or .fastq! Unrecognized file suffix for {0}! Abort!".format(file)

def make_sam_by_isoforms(isoform_fa_fq, read_stat_filename, flnc_fasta, output_dir='per_isoform_sam', nfl_fasta=None):
    fl_ids = set([r.id for r in SeqIO.parse(open(flnc_fasta), 'fasta')])
    nfl_ids = set([r.id for r in SeqIO.parse(open(nfl_fasta), 'fasta')]) if nfl_fasta is not None else None

    members = defaultdict(lambda: [])  # ex: PB.1.1 --> [list of fl members]

    print >> sys.stderr, "Reading read stat file", read_stat_filename
    for r in DictReader(open(read_stat_filename), delimiter='\t'):
        if r['pbid'] == 'NA': continue
        if r['id'] in fl_ids or (nfl_ids is not None and r['id'] in nfl_ids):
            members[r['pbid']].append(r['id'])

    print >> sys.stderr, "Lazy loading of FL fasta {0}...".format(flnc_fasta)
    fl_fd = SeqReaders.LazyFastaReader(flnc_fasta)
    if nfl_fasta is None:
        nfl_fd = None
    else:
        print >> sys.stderr, "Lazy loading of nFL fasta {0}...".format(nfl_fasta)
        nfl_fd = SeqReaders.LazyFastaReader(nfl_fasta)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    else:
        print >> sys.stderr, "WARNING: output directory {0} already exists. Overwriting!".format(output_dir)

    for r in SeqIO.parse(open(isoform_fa_fq), type_fa_or_fq(isoform_fa_fq)):
        pbid = r.id.split('|')[0]
        print >> sys.stderr, "Making SAM file for", pbid
        ref = open(os.path.join(output_dir, pbid + '_ref.fa'), 'w')
        inf = open(os.path.join(output_dir, pbid + '_in.fa'), 'w')
        ref.write(">{0}\n{1}\n".format(r.id, r.seq))
        ref.close()

        count = 0
        for x in members[pbid]:
            if x in fl_fd.keys():
                rec = fl_fd[x]
            else:
                rec = nfl_fd[x]

            ## in SA3 blasr will only allow output as SAM/BAM if the input records are <movie>/<start>_<end> format. so hack!
            if rec.id.endswith('_CCS'): rec.id = rec.id[:-4]
            inf.write(">{0}\n{1}\n".format(rec.id, rec.seq))
            count += 1
        inf.close()

        print pbid, count
        if count == 0:
            os.remove(ref.name)
            os.remove(inf.name)
            continue
        if sp_version == 'NEW':
            sp.blasr_for_quiver(inf.name, ref.name, os.path.join(output_dir, pbid + '.bam'),
                                bam = False, run_cmd = True, blasr_nproc = 12)
        else:
            sp.blasr_sam_for_quiver(inf.name, ref.name, os.path.join(output_dir, pbid + '.sam'))


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser("Make FL (+nFL) SAM files per final isoform. For debugging purpose only.")
    parser.add_argument("isoform_fa_fq", help="Final isoform fasta or fastq filename")
    parser.add_argument("read_stat", help="Read stat filename")
    parser.add_argument("--flnc_fa", default="isoseq_flnc.fasta", help="FLNC fasta filename (default: isoseq_flnc.fasta)")
    parser.add_argument("--nfl_fa", default=None, help="nFL fasta filename (ignored if not provided, default: OFF)")
    parser.add_argument("-o", "--output_dir", default="per_isoform_sam", help="Output directory (default: isoforms_per_sam/)")


    args = parser.parse_args()
    make_sam_by_isoforms(args.isoform_fa_fq, args.read_stat, args.flnc_fa, args.output_dir, args.nfl_fa)