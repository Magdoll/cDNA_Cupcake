"""
Experimemtal code for trimming primers & polyA tails from high error rate long reads
"""
import os, sys, pdb
from csv import DictWriter
from collections import namedtuple
from multiprocessing import Process
from Bio.Seq import Seq
from Bio import SeqIO
import parasail

ScoreTuple = namedtuple('ScoreTuple', ['score5', 'end5', 'score3', 'end3', 'endA'])

# for ONT using Clontech
#SEQ_5P = 'AAGCAGTGGTATCAACGCAGAGTACATGGGG'
#SEQ_3P_REV = 'GTATCAACGCAGAGTAC'

ISOSEQ_5P = 'GCAATGAAGTCGCAGGGTTGGG'
ISOSEQ_3P = 'GTACTCTGCGTTGATACCACTGCTT'

#SEQ_5P = 'GCAATGAAGTCGCAGGGTTGGGG'
#SEQ_5P = 'CAGGAAACAGCTATGACC'
#SEQ_3P_REV = 'AAGCAGTGGTATCAACGCAGAGTAC'
#SEQ_3P_REV = 'ACTGGCCGTCGTTTTAC'

MINSCORE_5P = 20
MINSCORE_3P = 20
MIN_A_LEN = 20

SCOREMAT = parasail.matrix_create("ACGT", 2, -5)

def trim5p3p_helper(r, seq_5p, seq_3p_rev):
    """
    Search for 5' and 3' in the first and last 100 bp window
    """
    s1 = str(r.seq[:100])
    s2 = str(r.reverse_complement().seq[:100])

    o1 = parasail.sg_qx_trace(s1, seq_5p, 3, 1, SCOREMAT)
    o2 = parasail.sg_qe_db_trace(s2, seq_3p_rev, 3, 1, SCOREMAT)
    lenA = None
    if o2.score >= MINSCORE_3P:
        lenA = trimA(s2[o2.end_query + 1:])

    if MIN_A_LEN == 0:
        end3 = len(r.seq) - o2.end_query - 1
        return ScoreTuple(score5=o1.score, end5=o1.end_query, score3=o2.score, end3=end3, endA=end3)
    elif lenA is not None:
        end3 = len(r.seq) - o2.end_query - 1
        endA = end3 - lenA + 1
        return ScoreTuple(score5=o1.score, end5=o1.end_query, score3=o2.score, end3=end3, endA=endA)
    else:
        end3 = len(r.seq) - o2.end_query - 1
        return ScoreTuple(score5=o1.score, end5=o1.end_query, score3=o2.score, end3=end3, endA=end3)

def trimA(rev_seq):
    if len(rev_seq) == 0:
        return None
    n_rev_seq = len(rev_seq)
    mismatch = 0
    i = 0
    while mismatch < 2 and i < n_rev_seq:
        if rev_seq[i]!='T':
            mismatch += 1
        i += 1
    i -= 1
    if i >= MIN_A_LEN:
        return i
    else:
        return None

def trim5p3p_multithreaded(fastq_filename, output_prefix, seq_5p, seq_3p_rev, chunks):
    # first figure out how many records there are and record positions
    num_lines = 0
    for line in open(fastq_filename, 'r'): num_lines += 1
    num_records = num_lines // 4
    chunk_size = (num_records//chunks) + (num_records%chunks>0)
    print("{0} has {1} records, {2} per chunk".format(fastq_filename, num_records, chunk_size))

    pools = []
    records = []
    count = 0
    i = 1
    for r in SeqIO.parse(open(fastq_filename), 'fastq'):
        count += 1
        records.append(r)
        if count >= chunk_size:
            p = Process(target=trim5p3p, args=(records, output_prefix+'.'+str(i), seq_5p, seq_3p_rev))
            p.start()
            print("Starting worker {0}...".format(i))
            pools.append(p)
            records = []
            count = 0
            i += 1
    p = Process(target=trim5p3p, args=(records, output_prefix + '.' + str(i), seq_5p, seq_3p_rev))
    p.start()
    print("Starting worker {0}...".format(i))
    pools.append(p)

    for p in pools:
        p.join()

    # now combine all the files
    f_FL = open(output_prefix+'.fl.fasta', 'w')
    f_FL_clips = open(output_prefix+'.fl.clips', 'w')
    f_nFL = open(output_prefix+'.nfl.fasta', 'w')
    f_csv = open(output_prefix+'.csv', 'w')
    for j in range(1, i+1):
        p = output_prefix + '.' + str(j)
        with open(p + '.fl.fasta') as h:
            f_FL.write(h.read())
            print("writing {0} into {1}...".format(h.name, f_FL.name))
        with open(p + '.fl.clips') as h:
            f_FL_clips.write(h.read())
            print("writing {0} into {1}...".format(h.name, f_FL_clips.name))
        with open(p + '.nfl.fasta') as h:
            f_nFL.write(h.read())
            print("writing {0} into {1}...".format(h.name, f_nFL.name))
        with open(p + '.csv') as h:
            f_csv.write(h.read())
            print("writing {0} into {1}...".format(h.name, f_csv.name))
        os.remove(p + '.fl.fasta')
        os.remove(p + '.fl.clips')
        os.remove(p + '.nfl.fasta')
        os.remove(p + '.csv')

    f_csv.close()
    f_FL.close()
    f_FL_clips.close()
    f_nFL.close()

def trim5p3p(records, output_prefix, seq_5p, seq_3p_rev):
    f_FL = open(output_prefix+'.fl.fasta', 'w')
    f_FL_clips = open(output_prefix+'.fl.clips', 'w')
    f_nFL = open(output_prefix+'.nfl.fasta', 'w')
    f_csv = open(output_prefix+'.csv', 'w')
    writer = DictWriter(f_csv, fieldnames=['id', 'end5', 'end3', 'endA', 'strand'])
    writer.writeheader()

    for r in records:
        r2 = r.reverse_complement()
        r2.id = r.id
        t1 = trim5p3p_helper(r, seq_5p, seq_3p_rev)
        t2 = trim5p3p_helper(r2, seq_5p, seq_3p_rev)

        is_fl_flag1 = t1.score5 >= MINSCORE_5P and t1.score3 >= MINSCORE_3P and (MIN_A_LEN == 0 or t1.endA!=t1.end3)
        is_fl_flag2 = t2.score5 >= MINSCORE_5P and t2.score3 >= MINSCORE_3P and (MIN_A_LEN == 0 or t2.endA!=t2.end3)

        if is_fl_flag1:
            if is_fl_flag2:
                if t1.score5+t1.score3 > t2.score5+t2.score3:
                    strand = '+'
                else:
                    strand = '-'
            else: # pick t1
                strand = '+'
        elif is_fl_flag2:
            strand = '-'
        else:
            strand = 'NA'

        info = {'id': r.id,
                'end5': 'NA',
                'end3': 'NA',
                'endA': 'NA',
                'strand': 'NA'}

        if strand == '+':
            info['strand'] = '+'
            info['end5'] = t1.end5
            info['end3'] = t1.end3
            info['endA'] = t1.endA
            f_FL.write(">{0}\n{1}\n".format(r.id, r.seq[t1.end5:t1.endA]))
            f_FL_clips.write(">{0}_5p strand:+ score:{1}\n{2}\n".format(r.id, t1.score5, r.seq[:t1.end5]))
            f_FL_clips.write(">{0}_3p strand:+ score:{1}\n{2}\n".format(r.id, t1.score3, r.seq[t1.endA:]))
        elif strand == '-':
            info['strand'] = '-'
            info['end5'] = t2.end5
            info['end3'] = t2.end3
            info['endA'] = t2.endA
            f_FL.write(">{0}\n{1}\n".format(r2.id, r2.seq[t2.end5:t2.endA]))
            f_FL_clips.write(">{0}_5p strand:- score:{1}\n{2}\n".format(r.id, t2.score5, r2.seq[:t2.end5]))
            f_FL_clips.write(">{0}_3p strand:- score:{1}\n{2}\n".format(r.id, t2.score3, r2.seq[t2.endA:]))
        else:
            # non-fL, but we still wanna trim away the stuff
            if t1.score5+t1.score3 > t2.score5+t2.score3:
                f_nFL.write(">{0} strand:+?\n{1}\n".format(r.id, r.seq[t1.end5:t1.endA]))
            else:
                f_nFL.write(">{0} strand:-?\n{1}\n".format(r2.id, r2.seq[t2.end5:t2.endA]))
        writer.writerow(info)
    f_csv.close()
    f_FL.close()
    f_FL_clips.close()
    f_nFL.close()


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("fastq_filename")
    parser.add_argument("output_prefix")
    parser.add_argument("-p", "--primer_fasta", default=None, help="Primer fasta file (if not given, use IsoSeq defaults)")
    parser.add_argument("-n", "--chunks", default=10, type=int, help="Number of chunks (CPUs) to use, default 10")

    args = parser.parse_args()

    if args.primer_fasta is None:
        seq_5p = ISOSEQ_5P
        seq_3p = ISOSEQ_3P
        print(f"Using Iso-Seq default 5' primer sequence: {seq_5p}")
        print(f"Using Iso-Seq default 3' primer sequence: {seq_3p}")
    else:
        reader = SeqIO.parse(open(args.primer_fasta), 'fasta')
        r = next(reader)
        if r.seqid!='5p':
            print("ERROR: the first entry in {0} should be >5p! Abort!".format(args.primer_fasta))
            sys.exit(-1)
        seq_5p = str(r.seq)
        r = next(reader)
        if r.seqid!='3p':
            print("ERROR: the second entry in {0} should be >3p! Abort!".format(args.primer_fasta))
            sys.exit(-1)
        seq_3p = str(r.seq)
        print(f"Reading in 5' primer sequence: {seq_5p}")
        print(f"Reading in 3' primer sequence: {seq_3p}")

    seq_3p_rev = str(Seq(seq_3p).reverse_complement())
    trim5p3p_multithreaded(args.fastq_filename, args.output_prefix, seq_5p, seq_3p_rev, args.chunks)
