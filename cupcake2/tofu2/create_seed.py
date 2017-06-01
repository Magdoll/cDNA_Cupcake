__author__ = 'lachesis'

import os, sys
from Bio import SeqIO
from cupcake.io.SeqReaders import LazyFastaReader
from cupcake2.io.FileIO import write_seqids_to_fasta

input = 'isoseq_flnc.fasta'

NUM_SEQS_PER_BATCH = 200000

d = LazyFastaReader(input)

lens = [(r.id, len(r.seq)) for r in SeqIO.parse(open(input), 'fasta')]
lens.sort(key=lambda x: x[1], reverse=True)

n = len(lens)

# start at 1% of the data
starting_seed_index = n * 1 / 100
good = [x[0] for x in lens[starting_seed_index:starting_seed_index+NUM_SEQS_PER_BATCH]]
write_seqids_to_fasta(good, 'seed0.fasta', d)

batch_index = 1
starting_index = starting_seed_index+NUM_SEQS_PER_BATCH
while starting_index < n:
    write_seqids_to_fasta([x[0] for x in lens[starting_index:starting_index+NUM_SEQS_PER_BATCH]], \
                          "batch{0}.fasta".format(batch_index), d)
    starting_index += NUM_SEQS_PER_BATCH
    batch_index += 1

write_seqids_to_fasta([x[0] for x in lens[:starting_seed_index]], "batch{0}.fasta".format(batch_index), d)
