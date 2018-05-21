
from csv import DictReader
from Bio import SeqIO
from ssw_wrap import Aligner

MIN_BARCODE_MATCH_LEN = 10  # have to at least see 10 of the 16 bp aligned
MIN_BARCODE_MATCH_SCORE = 20
MIN_BARCODE_SCORE_LEAD = 10

barcodes = ['atgacgcatcgtctga', 'gcagagtcatgtatag', 'gagtgctactctagta', 'catgtactgatacaca']
for i in xrange(4): barcodes[i] = barcodes[i].upper()
aligners = [Aligner(barcodes[i], match=2, mismatch=5, gap_open=3, gap_extend=1, report_secondary=False, report_cigar=True) for i in xrange(4)]


def main(ccs_fasta, flnc_fasta, primer_csv):
    good_flnc = []
    d = {}
    reader = DictReader(open(primer_csv),delimiter=',')
    for r in reader:
        zmw = r['id'][:r['id'].rfind('/')]
        d[zmw] = r

    flog = open(flnc_fasta+'.filtered.log', 'w')
    flog.write("flnc_id\tthreelen\tthreeseq\tscore0\tscore1\tscore2\tscore3\n")
    for r in SeqIO.parse(open(ccs_fasta), 'fasta'):
        zmw = r.id[:r.id.rfind('/')]
        if zmw not in d: continue # discarded short sequences
        rec = d[zmw]
        cands = [] # list of (barcode i, score, end-start)
        if rec['threeseen']=='1' and rec['fiveseen']=='1' and rec['polyAseen']=='1' and rec['chimera']=='0':
            # is FLNC
            end3 = int(rec['threeend'])
            flog.write("{0}\t{1}\t{2}".format(rec['id'], len(r.seq)-end3, r.seq[end3:].tostring()))
            for i,aligner in enumerate(aligners):
                o = aligner.align(r.seq[end3:].tostring())
                cands.append((i, o.score, o.ref_end-o.ref_begin))
                flog.write("\t{0}".format(o.score))
            flog.write("\n")
            cands.sort(key=lambda x: x[1], reverse=True)
            if cands[0][2] >= MIN_BARCODE_MATCH_LEN and cands[0][1] >= MIN_BARCODE_MATCH_SCORE and \
                    (len(cands)==1 or cands[0][1]-cands[1][1]>=MIN_BARCODE_SCORE_LEAD):
                good_flnc.append(rec['id'])
    flog.close()

    f = open(flnc_fasta + '.filtered.fasta', 'w')
    for r in SeqIO.parse(open(flnc_fasta),'fasta'):
        if r.id in good_flnc:
            f.write(">{0}\n{1}\n".format(r.description, r.seq))
    f.close()

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("ccs_fasta")
    parser.add_argument("flnc_fasta")
    parser.add_argument("primer_csv")

    args = parser.parse_args()
    main(args.ccs_fasta, args.flnc_fasta, args.primer_csv)