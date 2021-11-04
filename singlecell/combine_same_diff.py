import pysam
from csv import DictReader, DictWriter
from collections import defaultdict
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--same', '-s', help='Prefix for same, eg. output.same.deconcat.tagged')
parser.add_argument('--diff', '-d', help='Prefix for diff, eg. output.diff.deconcat.tagged')
parser.add_argument('--out', '-o', help='Prefix for output, eg. combined.flt')
args = parser.parse_args()

# note: remember to use the "old_id" 
# ex: m64182_210318_053335/0/ccs/1
in_diff = defaultdict(lambda: set()) # old zmw --> new zmw
in_same = defaultdict(lambda: set())

touse_diff = set()
touse_same = set()
# id      old_id  clip_len        extra   UMI     BC      BC_rev  BC_match        BC_top_rank     TSO
fields = ['id', 'old_id', 'clip_len', 'extra', 'UMI', 'BC', 'BC_rev', 'BC_match', 'BC_top_rank', 'TSO']
# f1 = open('combined.flt.csv', 'w')
f1 = open(args.out + '.csv', 'w')
writer = DictWriter(f1, fields, delimiter=',')
writer.writeheader()

# for r in DictReader(open('output.diff.deconcat.tagged.csv'),delimiter='\t'):
for r in DictReader(open(args.diff + '.csv'), delimiter='\t'):
    raw = r['old_id'].split('/')
    zmw = raw[0] + '/' + raw[1]
    in_diff[zmw].add(r['id'])
    touse_diff.add(r['id'])
    writer.writerow(r)

dup = 0    
# for r in DictReader(open('output.same.deconcat.tagged.csv'),delimiter='\t'):
for r in DictReader(open(args.same + '.csv'),delimiter='\t'):
    raw = r['old_id'].split('/')
    zmw = raw[0] + '/' + raw[1]
    if zmw in in_diff: dup += 1
    else:
        writer.writerow(r)
        in_same[zmw].add(r['id'])
        touse_same.add(r['id'])
    
f1.close()
print(f"{dup} zmws showed up in both diff/same...")


def get_zmw(seqid):
    raw = seqid.split('/')
    zmw = raw[0] + '/' + raw[1]
    return zmw

# reader = pysam.AlignmentFile(open('output.diff.deconcat.tagged.bam'),'rb',check_sq=False)
reader = pysam.AlignmentFile(open(args.diff + '.bam'),'rb',check_sq=False)
# f = pysam.AlignmentFile(open('combined.flt.bam', 'w'), 'wb', header=reader.header)
f = pysam.AlignmentFile(open(args.out + '.bam', 'w'), 'wb', header=reader.header)
for r in reader:
    if r.qname in touse_diff:
        print(r.qname)
        f.write(r)

# for r in pysam.AlignmentFile(open('output.same.deconcat.tagged.bam'),'rb',check_sq=False):
for r in pysam.AlignmentFile(open(args.same + '.bam'),'rb',check_sq=False):
    if r.qname in touse_same: f.write(r)

f.close()

