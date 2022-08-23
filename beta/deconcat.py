import os, sys, copy, subprocess
from csv import DictReader, DictWriter
from Bio import SeqIO
from Bio.Seq import Seq
import parasail
import pysam
from multiprocessing import Process
import pdb

SCOREMAT = parasail.matrix_create("ACGT", 2, -5)
MIN_SCORE = 60
MIN_LEN = 50  # minimum length post deconcatenation
MIN_PRE_LEN = 150 # minimum length pre-deconcatenation
PADDING_LEN = 3 # add a bit of buffer around the deconcatenated strings

global SEQ_R5_F5
global SEQ_F3_R3
global SEQ_F3_F5
global SEQ_R5_R3
global SEQ_METHOD

# for Jason's 10X USER 5' lib method
#SEQ_R5_F5='TCGGAAGAGCGTCGTGTAGATCGATCTACACGACGCTCTTCCGATCT'
#SEQ_F3_R3='GTACTCTGCGTTGATACCACTGCTTAGCGCTAAGCAGTGGTATCAACGCAGAGTAC'
#SEQ_R5_R3='TCGGAAGAGCGTCGTGTAGAGCGCTAAGCAGTGGTATCAACGCAGAGTAC '
#SEQ_F3_F5='GTACTCTGCGTTGATACCACTGCTTATCGATCTACACGACGCTCTTCCGATCT'

# for Jason's 10X USER 3' lib method
#SEQ_R5_F5='CCCATGTACTCTGCGTTGATACCACTGCTTATCGATAAGCAGTGGTATCAACGCAGAGTACATGGG'
#SEQ_R5_R3='CCCATGTACTCTGCGTTGATACCACTGCTTATCGATCTACACGACGCTCTTCCGATCT'
#SEQ_F3_R3='AGATCGGAAGAGCGTCGTGTAGAGCGCTCTACACGACGCTCTTCCGATCT'
#SEQ_F3_F5='AGATCGGAAGAGCGTCGTGTAGAGCGCTAAGCAGTGGTATCAACGCAGAGTACATGGG'

#SEQ_R5_F5='CCCCAACCCTGCGACTTCATTGCGCAATGAAGTCGCAGGGTTGGGG'
#SEQ_R5_R3='CCCCAACCCTGCGACTTCATTGCAAGCAGTGGTATCAACGCAGAGTAC'
#SEQ_F3_R3='GTACTCTGCGTTGATACCACTGCTTAAGCAGTGGTATCAACGCAGAGTAC'
#SEQ_F3_F5='GTACTCTGCGTTGATACCACTGCTTGCAATGAAGTCGCAGGGTTGGGG'

CSV_FIELDS = ['zmw','split','length','flag']



def deconcat_worker(input_bam, offset_start, offset_end, output_prefix, info, verbose=False):
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
        #print("rqname:",r.qname)
        zmw = r.qname[:r.qname.rfind('/')]
        start_flag = info[zmw]
        debug_flag = False
        #if r.qname.startswith('m64182_210416_060945/132915'): debug_flag = True
        it = deconcat_all(r.query, start_flag, start_pos=0, debug=debug_flag, verbose=verbose)
        i = 1
        for (s, e, flag, cur_seq) in it:
            if e-s < MIN_LEN:
                continue
            rec = {'zmw':zmw, 'split':i, 'length':e-s, 'flag':flag}
            writer.writerow(rec)
            #print(zmw,i,e-s,flag)
            d2 = copy.deepcopy(d)
            d2['name'] = d2['name'] + '/' + str(i)
            assert flag in ('F5', 'R3')
            d2['seq'] = d2['seq'][s:e]
            d2['qual'] = d2['qual'][s:e]
            if (SEQ_METHOD=='Jason-10X-5' and flag == 'R3') or \
                    (SEQ_METHOD=='Jason-10X-3' and flag == 'R3'):
                d2['seq'] = str(Seq(d2['seq']).reverse_complement())
                d2['qual'] = d2['qual'][::-1]
            x = pysam.AlignedSegment.from_dict(d2, r.header)
            f2.write(x)
            i += 1
    f1.close()
    f2.close()
    reader.close()



def deconcat_all(sequence, start_flag, start_pos, debug=False, verbose=False):
    cur_flag = start_flag
    cur_seq = sequence
    cur_offset = start_pos

    while len(cur_seq) > MIN_PRE_LEN:
        out = deconcat(cur_seq, cur_flag, debug, verbose)
        if out is None:
            if debug:
                pdb.set_trace()
            yield cur_offset, cur_offset+len(cur_seq), cur_flag, cur_seq
            break
        else:
            s, e, score, flag = out
            if debug: pdb.set_trace()
            # we need to see if there's an earlier match, with possibly just slightly worse score
            if s > MIN_PRE_LEN:
                alt_out = deconcat(cur_seq[:s], cur_flag, debug=debug, verbose=verbose)
            else:
                alt_out = None
            while alt_out is not None:
                alt_s, alt_e, alt_score, alt_flag = alt_out
                if alt_s > MIN_PRE_LEN:
                    earlier_alt_out = deconcat(cur_seq[:alt_s], cur_flag, debug=debug, verbose=verbose)
                else:
                    earlier_alt_out = None
                if earlier_alt_out is None: # we've reached as early as hit as we can, abort
                    break
                alt_out = earlier_alt_out
            if alt_out is not None:
                s, e, score, flag = alt_out

            s = min(len(cur_seq), s + PADDING_LEN)
            e = e - PADDING_LEN
            if cur_seq[:s] != sequence[cur_offset:(cur_offset+s)]:
                print(s)
                print(cur_seq[:s])
                print(sequence[cur_offset:(cur_offset+s)])
            assert cur_seq[:s] == sequence[cur_offset:(cur_offset+s)]
            yield cur_offset, cur_offset+s, cur_flag, cur_seq[:s]
            cur_seq = cur_seq[e:]
            cur_flag = flag
            cur_offset += e


def deconcat(sequence, prev, debug=False, verbose=False):
    if prev == 'R3':
        o1 = parasail.sg_qx_trace(sequence, SEQ_R5_F5, 3, 1, SCOREMAT)
        o2 = parasail.sg_qx_trace(sequence, SEQ_R5_R3, 3, 1, SCOREMAT)
        if o1.score >= MIN_SCORE and o1.score > o2.score:
            s1, e1 = o1.get_traceback().comp.find('|'), o1.end_query+1
            if verbose:
                print("DEBUG: F5", sequence[s1:e1])
            return s1, e1, o1.score, 'F5'
        elif o2.score >= MIN_SCORE:
            s2, e2 = o2.get_traceback().comp.find('|'), o2.end_query+1
            if verbose:
                print("DEBUG: R3", sequence[s2:e2])
            return s2, e2, o2.score, 'R3'
        else:
            if debug:
                pdb.set_trace()
            return None
    elif prev == 'F5':
        o1 = parasail.sg_qx_trace(sequence, SEQ_F3_R3, 3, 1, SCOREMAT)
        o2 = parasail.sg_qx_trace(sequence, SEQ_F3_F5, 3, 1, SCOREMAT)
        if o1.score >= MIN_SCORE and o1.score > o2.score:
            s1, e1 = o1.get_traceback().comp.find('|'), o1.end_query + 1
            if verbose:
                print("DEBUG: R3", sequence[s1:e1])
            return s1, e1, o1.score, 'R3'
        elif o2.score >= MIN_SCORE:
            s2, e2 = o2.get_traceback().comp.find('|'), o2.end_query+1
            if verbose:
                print("DEBUG: F5", sequence[s2:e2])
            return s2, e2, o2.score, 'F5'
        else:
            if debug:
                pdb.set_trace()
            return None
    else:
        raise ValueError("Expected previous primer to be F5 or R3. Saw {0} instead. Abort!".format(prev))


def main(input_prefix, output_prefix, cpus, verbose=False):
    info = {}
    for r in SeqIO.parse(open(input_prefix + '.lima.clips'), 'fasta'):
        zmw = r.id[:r.id.rfind('/')]
        e = int(r.id.split('/')[2].split('_')[1])
        if e < 100:
            info[zmw] = 'F5' if r.description.split('bc:')[-1] == '0' else 'R3'
    print("Finished reading lima clips file.", file=sys.stdout)

    num_records = len(info)
    input_bam = input_prefix + '.bam'
    onames = []
    if cpus == 1:
        oname = output_prefix
        deconcat_worker(input_bam, 0, num_records, oname, info, verbose)
    else:
        chunk_size = (num_records // cpus) + (num_records % cpus)
        offset_start = 0
        pools = []
        while offset_start < num_records:
            oname = output_prefix+'.'+str(offset_start)
            p = Process(target=deconcat_worker, args=(input_bam, offset_start, offset_start+chunk_size, oname, info, verbose))
            p.start()
            print("Launching deconcat worker for records {0}-{1}...".format(offset_start, offset_start+chunk_size))
            offset_start += chunk_size
            pools.append(p)
            onames.append(oname)
        for p in pools:
            print("Waiting for {0} to finish...".format(p.name))
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
    CMD = "bamtools merge -in " + " -in ".join(bams) + " -out " + output_prefix + '.bam'
    subprocess.check_call(CMD, shell=True)

    for oname in onames:
        os.remove(oname+'.bam')
        os.remove(oname+'.csv')
    print("Output written to: {0}.bam, {0}.csv".format(output_prefix))


if __name__ == "__main__":
    global SEQ_R5_F5
    global SEQ_F3_R3
    global SEQ_F3_F5
    global SEQ_R5_R3
    global SEQ_METHOD

    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("input_prefix")
    parser.add_argument("output_prefix")
    parser.add_argument("-n", "--cpus", type=int, default=10, help="Number of CPUs")
    parser.add_argument("-m", "--method", choices=['Jason-10X-5', 'Jason-10X-3'])
    parser.add_argument("--verbose", default=False, action="store_true")

    args = parser.parse_args()

    if args.method == 'Jason-10X-5':
        SEQ_R5_F5='TCGGAAGAGCGTCGTGTAGATCGATCTACACGACGCTCTTCCGATCT'
        SEQ_F3_R3='GTACTCTGCGTTGATACCACTGCTTAGCGCTAAGCAGTGGTATCAACGCAGAGTAC'
        SEQ_R5_R3='TCGGAAGAGCGTCGTGTAGAGCGCTAAGCAGTGGTATCAACGCAGAGTAC '
        SEQ_F3_F5='GTACTCTGCGTTGATACCACTGCTTATCGATCTACACGACGCTCTTCCGATCT'
        SEQ_METHOD = args.method
    elif args.method == 'Jason-10X-3':
        SEQ_R5_F5='CCCATGTACTCTGCGTTGATACCACTGCTTATCGATAAGCAGTGGTATCAACGCAGAGTACATGGG'
        SEQ_R5_R3='CCCATGTACTCTGCGTTGATACCACTGCTTATCGATCTACACGACGCTCTTCCGATCT'
        SEQ_F3_R3='AGATCGGAAGAGCGTCGTGTAGAGCGCTCTACACGACGCTCTTCCGATCT'
        SEQ_F3_F5='AGATCGGAAGAGCGTCGTGTAGAGCGCTAAGCAGTGGTATCAACGCAGAGTACATGGG'
        SEQ_METHOD = args.method

    main(args.input_prefix, args.output_prefix, args.cpus, args.verbose)
