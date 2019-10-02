__author__ = 'etseng@pacb.com'

"""
Given a list of directories (ex: preCluster_out/0, preCluster_out/1...)
that have already run through IceIterative2,

1. collect output/final.consensus.fasta - modify the cids
2. collect output/final.pickle - modify the cids

Modify the cids from "c0" --> "pC<bin>_c0"
"""

import os, sys, glob
import fileinput
from Bio import SeqIO
from pickle import *
from cupcake.io.SeqReaders import LazyFastaReader


def collect_ice2_dirs(dirs_to_collect, output_fasta, output_pickle, new_cid_format="b{bin}_c{cid}"):

    outf = open(output_fasta, 'w')
    combined_uc = {}
    combined_refs = {}


    for d1 in dirs_to_collect:
        bin = os.path.basename(d1)
        fasta_file = os.path.join(d1, 'output', 'final.consensus.fasta')
        pickle_file = os.path.join(d1, 'output', 'final.pickle')
        a = load(open(pickle_file))
        uc = a['uc']
        refs = a['refs']
        for cid in uc:
            newcid = new_cid_format.format(bin=bin, cid=cid)
            combined_uc[newcid] = uc[cid]
            # also need to re-write the refs
            # in case where the refs no longer exist (user deletion, just SKIP)
            new_ref = os.path.join(os.path.dirname(refs[cid]), 'collected_ref.fasta')
            if os.path.exists(refs[cid]):
                r  = next(SeqIO.parse(open(refs[cid]), 'fasta'))
                with open(new_ref, 'w') as f:
                    f.write(">{newcid}\n{seq}\n".format(newcid=newcid, seq=r.seq))
            else:
                print("WARNING: ref file {0} no longer exists. new ref files not written.".format(refs[cid]), file=sys.stderr)
            combined_refs[newcid] = new_ref

        for line in fileinput.input(fasta_file):
            if line.startswith('>'):
                cid, rest = line.split('/', 1)
                # assert cid.startswith('>c')
                cid = int(cid[2:])
                outf.write(">{newcid}/{rest}".format(\
                    newcid=new_cid_format.format(bin=os.path.basename(d1), cid=cid), rest=rest))
            else:
                outf.write(line)

    outf.close()
    with open(output_pickle, 'w') as f:
        dump({'uc': combined_uc, 'refs': combined_refs}, f)

    return combined_uc, combined_refs


# def collect_ice2_dirs(dirs_to_collect, output_fasta, output_pickle, new_cid_format="b{bin}_c{cid}"):
#
#     outf = open(output_fasta, 'w')
#     combined_uc = {}
#     combined_refs = {}
#
#     for d1 in dirs_to_collect:
#         bin = os.path.basename(d1)
#         fasta_file = os.path.join(d1, 'output', 'final.consensus.fasta')
#         pickle_file = os.path.join(d1, 'output', 'final.pickle')
#         a = load(open(pickle_file))
#         uc = a['uc']
#         refs = a['refs']
#         for cid in uc:
#             newcid = new_cid_format.format(bin=bin, cid=cid)
#             combined_uc[newcid] = uc[cid]
#             # also need to re-write the refs
#             new_ref = os.path.join(os.path.dirname(refs[cid]), 'collected_ref.fasta')
#             r  = SeqIO.parse(open(refs[cid]), 'fasta').next()
#             with open(new_ref, 'w') as f:
#                 f.write(">{newcid}\n{seq}\n".format(newcid=newcid, seq=r.seq))
#             combined_refs[newcid] = new_ref
#         for r in SeqIO.parse(open(fasta_file), 'fasta'):
#             cid, rest = r.id.split('/', 1)
#             assert cid.startswith('c')
#             cid = int(cid[1:])
#             outf.write(">{newcid}/{rest}\n{seq}\n".format(\
#                 newcid=new_cid_format.format(bin=os.path.basename(d1), cid=cid), rest=rest, seq=r.seq))
#     outf.close()
#     with open(output_pickle, 'w') as f:
#         dump({'uc': combined_uc, 'refs': combined_refs}, f)
#
#     return combined_uc, combined_refs


def chunk_collected_fasta_pickle(combined_fasta, combined_uc, combined_refs, num_chunks, chunk_prefix):
    d = LazyFastaReader(combined_fasta, seqid_extraction=lambda x: x.split('/')[0])
    # the seqids in combined_fasta are b0_c0/2/1776, need to make b0_c0 also key

    uc = combined_uc
    refs = combined_refs

    keys = list(uc.keys())
    keys.sort()
    n = len(keys) / num_chunks + 1
    for i in range(num_chunks):
        _from = i * n
        _to = min(len(keys), (i+1) * n)
        with open("{0}.chunk{1}.consensus.fasta".format(chunk_prefix, i), 'w') as f:
            for seqid in keys[_from:_to]:
                r = d[seqid]
                f.write(">{0}\n{1}\n".format(r.id, r.seq))
        with open("{0}.chunk{1}.pickle".format(chunk_prefix, i), 'w') as f:
            dump({'uc': dict((k, uc[k]) for k in keys[_from:_to]),
                  'refs': dict((k, refs[k]) for k in keys[_from:_to])}, f)


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser("Collect ICE results (final.consensus.fasta and final.pickle)")
    parser.add_argument("ice_dirs_file",  help="File containing list of Ice2 directories")
    parser.add_argument("output_prefix", help="Output prefix (ex: final_collected)")
    parser.add_argument("--num_chunks", default=10, type=int, help="Number of chunks (default: 10)")

    args = parser.parse_args()
    dirs_to_collect = [line.strip() for line in open(args.ice_dirs_file)]

    output_fasta = args.output_prefix + '.consensus.fasta'
    output_pickle = args.output_prefix + '.pickle'

    combined_uc, combined_refs = collect_ice2_dirs(dirs_to_collect, output_fasta, output_pickle)

    chunk_collected_fasta_pickle(combined_fasta=output_fasta,
                                 combined_uc=combined_uc,
                                 combined_refs=combined_refs,
                                 num_chunks=args.num_chunks,
                                 chunk_prefix=args.output_prefix)