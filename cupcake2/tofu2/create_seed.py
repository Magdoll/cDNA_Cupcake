import os, sys
from Bio import SeqIO
from cupcake.io.SeqReaders import LazyFastaReader
from cupcake2.io.FileIO import write_seqids_to_fasta

def create_seed_n_batch_files(input='isoseq_flnc.fasta', fasta_d=None, seed_filename='seed0.fasta', batch_pre='batch', num_seqs_per_batch=100000):
    if fasta_d is None:
        fasta_d = LazyFastaReader(input)

    batch_files = []

    lens = [(r.id, len(r.seq)) for r in SeqIO.parse(open(input), 'fasta')]
    lens.sort(key=lambda x: x[1], reverse=True)

    n = len(lens)

    # start at 1% of the data
    starting_seed_index = n * 1 / 100
    good = [x[0] for x in lens[starting_seed_index:starting_seed_index+num_seqs_per_batch]]
    write_seqids_to_fasta(good, seed_filename, fasta_d)

    batch_index = 1
    starting_index = starting_seed_index+num_seqs_per_batch
    while starting_index < n:
        write_seqids_to_fasta([x[0] for x in lens[starting_index:starting_index+num_seqs_per_batch]], \
                              "{0}{1}.fasta".format(batch_pre, batch_index), fasta_d)
        starting_index += num_seqs_per_batch
        batch_index += 1
        batch_files.append("{0}{1}.fasta".format(batch_pre, batch_index))

    write_seqids_to_fasta([x[0] for x in lens[:starting_seed_index]], "{0}{1}.fasta".format(batch_pre, batch_index), fasta_d)
    return batch_index+1

if __name__ == "__main__":
    create_seed_n_batch_files()
