__author__ = 'etseng@pacb.com'

import random
from cupcake.io.SeqReaders import LazyFastaReader

def write_select_seqs_to_fasta(fasta_filename, seqids, output_filename, mode='w'):
    d = LazyFastaReader('isoseq_flnc.fasta')
    with open(output_filename, mode) as f:
        r = d[x]
        f.write(">{0}\n{1}\n".format(r.id, r.seq))


def write_preClusterSet_to_fasta(pCS, output_filename, fasta_d):
    """
    Write to fasta:
    ID -- cid | selected representative seqid for this cid
    Seq --- sequence of the selected representative

    Currently, the rep is randomly chosen.
    """
    with open(output_filename, 'w') as f:
        for cid in pCS.S:
            r = fasta_d[random.choice(pCS.S[cid].members)]
            f.write(">{0}\n{1}\n".format(r.id, r.seq))


def write_seqids_to_fasta(seqids, output_filename, fasta_d):
    """
    Write to fasta:
    ID --- the sequence id
    Seq -- the sequence
    """
    with open(output_filename, 'w') as f:
        for seqid in seqids:
            r = fasta_d[seqid]
            f.write(">{0}\n{1}\n".format(r.id, r.seq))
