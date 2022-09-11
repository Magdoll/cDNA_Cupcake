#!/usr/bin/env python

__version__ = '1.0'

import os, sys
import random
from collections import namedtuple, defaultdict
from Bio import SeqIO

simType = ['sub', 'ins', 'del', 'match']
simTypeSize = 4

def throwdice(profile):
    dice = random.random()
    for i in range(simTypeSize):
        if dice < profile[i]:
            return simType[i]
        
def sim_start(ntimes, profile):
    start = defaultdict(lambda: 0)
    acc = defaultdict(lambda: 0)
    for i in range(ntimes):
        curpos = 0
        while True:
            type = throwdice(profile)
            acc[type] += 1
            if type=='match': 
                start[curpos] += 1
                break
            elif type=='ins':
                # insertion occurred, with 1/4 chance it'll match
                if random.random() < .25:
                    start[curpos] += 1
                    break
            # if is sub or del, just advance cursor forward
            curpos += 1
    return start, acc

def sim_seq(seq, profile):
    """
    :param seq: sequence to simulate from
    :param profile: accumulative prob vector for simType, ex: [0.01, 0.05, 0.07, 1.] means 1% sub, 4% ins, 2% del
    :return: simulated sequence, qv string ('!' for err, ']' for no err, currently don't handle deletion)
    """
    nucl = set(['A','T','C','G'])
    sim = ''
    qv = ''  # qv string,
    for i, s in enumerate(seq):
        while True:
            type = throwdice(profile)
            if type=='match': 
                sim += s
                qv += ']'
                break
            elif type=='ins':
                # insertion occurred, with 1/4 chance it'll match
                choice = random.sample(nucl,1)[0]
                sim += choice
                qv += '!'
            elif type=='sub': # anything but the right one
                choice = random.sample(nucl.difference([s]),1)[0]
                sim += choice
                qv += '!'
                break
            elif type=='del': # skip over this
                break
            else: raise KeyError("Invalid type {0}".format(type))
        
    return sim, qv
            
if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser("Simple error simulation")
    parser.add_argument("fasta_filename")
    parser.add_argument("output_prefix")
    parser.add_argument("--copy", type=int, default=1, help="Number of copies to simulate per input sequence (default: 1)")
    parser.add_argument("--ins", "-i", type=float, default=0, help="Insert error rate [0-1] (default: 0)")
    parser.add_argument("--dele", "-d", type=float, default=0, help="Deletion error rate [0-1] (default: 0)")
    parser.add_argument("--sub", "-s", type=float, default=0, help="Substitution error rate [0-1] (default: 0)")

    args = parser.parse_args()

    if args.sub < 0 or args.sub > 1: 
        print("Substitution error must be between 0-1!", file=sys.stderr)
        sys.exit(-1)
    if args.ins < 0 or args.ins > 1:
        print("Insertion error must be between 0-1!", file=sys.stderr)
        sys.exit(-1)
    if args.dele < 0 or args.dele > 1:
        print("Deletion error must be between 0-1!", file=sys.stderr)
        sys.exit(-1)

    if args.sub + args.ins + args.dele > 1:
        print("Total sub+ins+del error cannot exceed 1!", file=sys.stderr)
        sys.exit(-1)


    profile = [args.sub, args.sub+args.ins, args.ins+args.dele, 1.]

    fasta_filename = args.fasta_filename
    idpre = args.output_prefix

    ith = 0
    for r in SeqIO.parse(open(fasta_filename), 'fasta'):
        for j in range(args.copy):
            ith += 1
            print((">{0}_{1}_{2}\n{3}".format(idpre, ith, r.id[:r.id.find('|')], sim_seq(r.seq.tostring(), profile))))
