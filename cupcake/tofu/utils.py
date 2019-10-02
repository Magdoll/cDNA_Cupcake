__author__ = 'etseng@pacb.com'

from Bio import SeqIO

def check_ids_unique(fa_or_fq_filename, is_fq=False):
    """
    Confirm that a FASTA/FASTQ file has all unique IDs
    (used probably by collapse or fusion finding script)
    """
    seen = set()
    for r in SeqIO.parse(open(fa_or_fq_filename), 'fastq' if is_fq else 'fasta'):
        if r.id in seen:
            raise Exception("Duplicate id {0} detected. Abort!".format(r.id))
        seen.add(r.id)
