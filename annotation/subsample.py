#!/usr/bin/env python
import os, sys, random
import math
from csv import DictReader
from collections import defaultdict
from multiprocessing import Pool

def get_counts(count_filename, min_fl_count=2, key='id', min_len=None, max_len=None):
    total = 0
    count_d = defaultdict(lambda: 0)
    for r in DictReader(open(count_filename), delimiter='\t'):
        _len = int(r['length'])
        if min_len is not None and _len < min_len: continue
        if max_len is not None and _len > max_len: continue
        c = int(r['fl_count'])
        if c >= min_fl_count:
            count_d[r[key]] += c
            total += c
    
    counts = []
    for k,v in count_d.items():
        counts += [k]*v
    
    return total, counts

def sampleBySize(s, iteration, counts):
    tmp = []
    for i in range(iteration):
        tally = defaultdict(lambda: 0)
        for k in random.sample(counts, s):
            tally[k] += 1
        tmp.append(len(tally)) 
    _mean = sum(tmp)*1./len(tmp)
    _std = math.sqrt(sum((x-_mean)**2 for x in tmp)*1./len(tmp))
    return((s, min(tmp), max(tmp), _mean, _std))

def subsample_parallel(total, counts, iteration=100, min_fl_count=2, step=10**4, ncore=1):
    sizes = list(range(0, total+1, step))
    print("#min fl count:", min_fl_count)
    print("size", "min", "max", "mean", "sd")
    
    with Pool(ncore) as p:
        return p.starmap(sampleBySize, [(s, iteration, counts) for s in sizes])

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("count_filename")
    parser.add_argument("--by", default='id', help="Unique specifier name(default: id)")
    parser.add_argument("--iterations", type=int, default=100, help="Number of iterations (default: 100)")
    parser.add_argument("--range", default=None, help="Length range (ex: (1000,2000), default None)")
    parser.add_argument("--min_fl_count", default=1, type=int, help="Minimum FL count (default: 1)")
    parser.add_argument("--step", default=10000, type=int, help="Step size (default: 10000)")
    parser.add_argument("--ncores", default=1, type=int, help="Number of threads to use for parallelization")
    args = parser.parse_args()

    min_len, max_len = None, None
    if args.range is not None:
        min_len, max_len = eval(args.range)
        assert 0 <= min_len < max_len

    total, counts = get_counts(args.count_filename, args.min_fl_count, args.by, min_len, max_len)
    subsample_parallel(total, counts, args.iterations, args.min_fl_count, args.step, args.ncores)


