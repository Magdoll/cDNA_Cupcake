#!/usr/bin/env python
import os, sys, random
import math
from csv import DictReader
from collections import defaultdict

def get_counts(count_filename, min_fl_count=2, key='id', min_len=None, max_len=None):
    total = 0
    count_d = defaultdict(lambda: 0)
    for r in DictReader(open(count_filename), delimiter='\t'):
        _len = int(r['length'])
        if min_len is not None and _len < min_len: continue
        if max_len is not None and _len > max_len: continue
        c = int(r['fl_count'])
        if c >= min_fl_count:
            count_d[r[key]+'---'+r['category']] += c
            total += c
    
    counts = []
    for k,v in count_d.items():
        counts += [k]*v
    
    return total, counts

from collections import defaultdict
def subsample(total, counts, iter=100, min_fl_count=2, step=10**4):
    sizes = list(range(100, total+1, step))
    print("min fl count:", min_fl_count)
    print("size", "category", "min", "max", "mean", "sd")
    for s in sizes:
        tmp = defaultdict(lambda: []) # category --> N iterations
        for i in range(iter):
            tally = defaultdict(lambda: 0)
            uniq_id_count = defaultdict(lambda: set()) # category -> unique ids
            for k in random.sample(counts, s):
                tally[k] += 1
            for k in tally:
                _id, _cat = k.split('---')
                uniq_id_count[_cat].add(_id) #if tally[k] >= min_fl_count: uniq_id_count[_cat].add(_id)
            for _cat in uniq_id_count:
                tmp[_cat].append(len(uniq_id_count[_cat]))

        for _cat in tmp:
            _mean = sum(tmp[_cat])*1./len(tmp[_cat])
            _std = math.sqrt(sum((x-_mean)**2 for x in tmp[_cat])*1./len(tmp[_cat]))
            print(s, _cat, min(tmp[_cat]), max(tmp[_cat]), _mean, _std)

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("count_filename")
    parser.add_argument("--by", default='id', help="Unique specifier name(default: id)")
    parser.add_argument("--iterations", type=int, default=100, help="Number of iterations (default: 100)")
    parser.add_argument("--range", default=None, help="Length range (ex: (1000,2000), default None)")
    parser.add_argument("--min_fl_count", default=1, type=int, help="Minimum FL count (default: 1)")
    parser.add_argument("--step", default=10000, type=int, help="Step size (default: 10000)")
    args = parser.parse_args()

    min_len, max_len = None, None
    if args.range is not None:
        min_len, max_len = eval(args.range)
        assert 0 <= min_len < max_len

    total, counts = get_counts(args.count_filename, args.min_fl_count, args.by, min_len, max_len)
    subsample(total, counts, args.iterations, args.min_fl_count, args.step)


