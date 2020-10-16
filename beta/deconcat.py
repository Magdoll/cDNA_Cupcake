import os, sys
import copy
from csv import DictReader, DictWriter
from Bio import SeqIO
from Bio.Seq import Seq
import parasail
import pysam
from multiprocessing import Process

SCOREMAT = parasail.matrix_create("ACGT", 2, -5)
MIN_SCORE = 80
MIN_LEN = 50

SEQ_R5_F5='CCCCAACCCTGCGACTTCATTGCGCAATGAAGTCGCAGGGTTGGGG'
SEQ_R5_R3='CCCCAACCCTGCGACTTCATTGCAAGCAGTGGTATCAACGCAGAGTAC'
SEQ_F3_R3='GTACTCTGCGTTGATACCACTGCTTAAGCAGTGGTATCAACGCAGAGTAC'
SEQ_F3_F5='GTACTCTGCGTTGATACCACTGCTTGCAATGAAGTCGCAGGGTTGGGG'

CSV_FIELDS = ['zmw','split','length','flag']

def deconcat_worker(input_bam, offset_start, offset_end, output_prefix, info):
    reader = pysam.AlignmentFile(input_bam, 'rb', check_sq=False)
    f1 = open(output_prefix + '.csv', 'w')
    writer = DictWriter(f1, CSV_FIELDS, delimiter=',')
    writer.writeheader()
    f2 = pysam.AlignmentFile(output_prefix + '.bam', 'wb', header=reader.header)
    counter = -1
    for r in reader:
        counter += 1
        if counter < offset_start: continue
        if counter >= offset_end: break

        d = r.to_dict()
        zmw = r.qname[:r.qname.rfind('/')]
        start_flag = info[zmw]
        it = deconcat_all(r.query, start_flag, start_pos=0)
        i = 1
        for (s, e, flag, cur_seq) in it:
            if e-s < MIN_SCORE:
                continue
            rec = {'zmw':zmw, 'split':i, 'length':e-s, 'flag':flag}
            writer.writerow(rec)
            #print(zmw,i,e-s,flag)
            d2 = copy.deepcopy(d)
            d2['name'] = d2['name'] + '/' + str(i)
            assert flag in ('F5', 'R3')
            d2['seq'] = d2['seq'][s:e]
            d2['qual'] = d2['qual'][s:e]
            if flag == 'R3':
                d2['seq'] = str(Seq(d2['seq']).reverse_complement())
                d2['qual'] = d2['qual'][::-1]

            x = pysam.AlignedSegment.from_dict(d2, r.header)
            f2.write(x)
            i += 1
    f1.close()
    f2.close()


def deconcat_all(sequence, start_flag, start_pos):
    cur_flag = start_flag
    cur_seq = sequence
    cur_pos = start_pos

    while len(cur_seq) > MIN_LEN:
        out = deconcat(cur_seq, cur_flag)
        if out is None:
            yield cur_pos, cur_pos+len(cur_seq), cur_flag, cur_seq
            break
        else:
            s, e, score, flag = out
            yield cur_pos, cur_pos+s, cur_flag, cur_seq[:s]
            cur_seq = cur_seq[e:]
            cur_flag = flag
            cur_pos = e


def deconcat(sequence, prev):
    if prev == 'R3':
        o1 = parasail.sg_qx_trace(sequence, SEQ_R5_F5, 3, 1, SCOREMAT)
        o2 = parasail.sg_qx_trace(sequence, SEQ_R5_R3, 3, 1, SCOREMAT)
        if o1.score >= MIN_SCORE and o1.score > o2.score:
            return o1.get_traceback().comp.find('|'), o1.end_query, o1.score, 'F5'
        elif o2.score >= MIN_SCORE:
            return o2.get_traceback().comp.find('|'), o2.end_query, o2.score, 'R3'
        else:
            return None
    elif prev == 'F5':
        o1 = parasail.sg_qx_trace(sequence, SEQ_F3_R3, 3, 1, SCOREMAT)
        o2 = parasail.sg_qx_trace(sequence, SEQ_F3_F5, 3, 1, SCOREMAT)
        if o1.score >= MIN_SCORE and o1.score > o2.score:
            return o1.get_traceback().comp.find('|'), o1.end_query, o1.score, 'R3'
        elif o2.score >= MIN_SCORE:
            return o2.get_traceback().comp.find('|'), o2.end_query, o2.score, 'F5'
        else:
            return None
    else:
        raise ValueError("Expected previous primer to be F5 or R3. Saw {0} instead. Abort!".format(prev))


def main(input_prefix, output_prefix, cpus):
    info = {}
    for r in SeqIO.parse(open(input_prefix + '.lima.clips'), 'fasta'):
        zmw = r.id[:r.id.rfind('/')]
        e = int(r.id.split('/')[2].split('_')[1])
        if e < 100:
            info[zmw] = 'F5' if r.description.split('bc:')[-1] == '0' else 'R3'
    print("Finished reading lima clips file.", file=sys.stdout)


    num_records = len(info)
    chunk_size = (num_records // cpus) + (num_records % cpus)

    offset_start = 0
    input_bam = input_prefix + '.bam'
    pools = []
    onames = []
    while offset_start <= num_records:
        oname = output_prefix+'.'+str(offset_start)
        p = Process(target=deconcat_worker, args=(input_bam, offset_start, offset_start+chunk_size, oname, info,))
        p.start()
        print("Launching deconcat worker for records {0}-{1}...".format(offset_start, offset_start+chunk_size))
        offset_start += chunk_size
        pools.append(p)
        onames.append(oname)

    for p in pools:
        p.join()

    print("All deconcat workers done. Collecting results.")
    f_csv = open(output_prefix + '.csv', 'w')
    writer = DictWriter(f_csv, CSV_FIELDS, delimiter=',')
    writer.writeheader()
    bams = []
    for oname in onames:
        bams.append(oname+'.bam')
        for r in DictReader(open(oname + '.csv'), delimiter=','):
            writer.writerow(r)
    f_csv.close()
    print("Merging bam files...")
    reader = pysam.AlignmentFile(bams[0], 'rb', check_sq=False)
    f = pysam.AlignmentFile(output_prefix + '.bam', 'wb', header=reader.header)
    for bam in bams:
        for r in pysam.AlignmentFile(bam, 'rb', check_sq=False):
            x = pysam.AlignedSegment.from_dict(r.to_dict(), r.header)
            f.write(x)
    f.close()

    #pysam.merge(output_prefix+'.bam', *bams)

    for oname in onames:
        os.remove(oname+'.bam')
        os.remove(oname+'.csv')
    print("Output written to: {0}.bam, {0}.csv".format(output_prefix))


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("input_prefix")
    parser.add_argument("output_prefix")
    parser.add_argument("-n", "--cpus", type=int, default=10, help="Number of CPUs")

    args = parser.parse_args()
    main(args.input_prefix, args.output_prefix, args.cpus)
