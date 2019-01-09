from Bio import SeqIO
from Bio.Seq import Seq
from csv import DictReader, DictWriter
import pysam
from pbcore.io import BamIO


def find_Aend(seq, min_a_len=8):
    """
    Given a sequence, find the likely beginning and end of polyA tail
    """
    Aseq = 'A'*min_a_len
    # will search for only the last 200 bp of the sequence
    x = seq[-200:]
    j = x.rfind(Aseq)

    if j >= 0:
        # now find the beginning
        end = len(seq)-200+j+min_a_len
        start = end-1
        while start >= 0 and seq[start]=='A':
            start -= 1
        return start, end
    else:
        return -1, -1

def find_Gstart(seq, min_g_len=3):
    """
    For Brendan's UMI design
    5' NEB --- (NNNNHHHH) --- GGG(G) --- transcript --- (A)n --- 3' NEB

    the sequence should already have NEB primers removed.
    detect the beginning and end of G's and return it
    """
    len_seq = len(seq)
    i = seq.find('G'*min_g_len)
    if i > 0:
        # now find the end
        j = i + 3
        while j < len_seq and seq[j]=='G':
            j += 1
        return i, j
    else:
        return -1, -1



def clip_out(bam_filename, umi_len, bc_len, output_prefix, UMI_type, shortread_bc={}):
    """
    :param bam_filename: BAM of post-LIMA (primer-trimmed) CCS sequences
    :param UMI_type: either 'A3' or 'G5'
    :param shortread_bc: a dict of barcode -> "Y|N" for top-ranked. If given, came from short read data.
    """
    assert UMI_type in ('A3', 'G5')

    umi_bc_len = umi_len + bc_len

    FIELDS = ['id', 'clip_len', 'extra', 'UMI', 'BC', 'BC_rev', 'BC_match', 'BC_top_rank']

    f1 = open(output_prefix + '.trimmed.csv', 'w')
    writer1 = DictWriter(f1, FIELDS, delimiter='\t')
    writer1.writeheader()

    reader = BamIO.IndexedBamReader(bam_filename)
    f2 = pysam.AlignmentFile(output_prefix+'.trimmed.bam', 'wb', template=reader.peer)

    for r in reader:
        d = r.peer.to_dict()

        if not r.isForwardStrand: # need to revcomp the seq and qual
            d['seq'] = str(Seq(d['seq']).reverse_complement())
            d['qual'] = d['qual'][::-1]
            new_tags = []
            for tag in d['tags']:
                if tag.startswith('dq:i:') or tag.startswith('iq:i:') or tag.startswith('sq:i:'):
                    tag = tag[:5] + tag[::-1][:-5]
                new_tags.append(tag)
            d['tags'] = new_tags
            d['flag'] = '4'   # convert it back to not being rev complemented


        if UMI_type == 'A3':
            A_start, A_end = find_Aend(d['seq'])
            if A_end > 0:
                seq2 = d['seq'][A_end:]  # should be just UMI + BC, unless UMI started with 'A's

                diff = len(seq2) - umi_bc_len
                if diff < 0: # UMI may have started with 'A's
                    seq2 = d['seq'][A_end+diff:]

                seq_extra = 'NA'
                if diff > 0: seq_extra = seq2[:diff]

                seq_bc = seq2[-bc_len:]
                seq_umi = seq2[-(bc_len+umi_len):-bc_len]


                # reverse complement BC because it's always listed in rev comp in short read data
                seq_bc_rev = str(Seq(seq_bc).reverse_complement())

                match = 'Y' if seq_bc_rev in shortread_bc else 'N'
                match_top = 'Y' if (match=='Y' and shortread_bc[seq_bc_rev]=='Y') else 'N'

                rec = {'id': r.peer.qname,
                       'clip_len': len(seq2),
                       'extra': seq_extra,
                       'UMI': seq_umi,
                       'BC': seq_bc,
                       'BC_rev': seq_bc_rev,
                       'BC_match': match,
                       'BC_top_rank': match_top}
                writer1.writerow(rec)


                # subset the sequence to include only the polyAs
                d['seq'] = d['seq'][:A_end]
                d['qual'] = d['qual'][:A_end]
                assert len(d['seq'])==len(d['qual'])
                new_tags = []
                for tag in d['tags']:
                    if tag.startswith('zs:B'): # defunct CCS tag, don't use
                        pass
                    elif tag.startswith('dq:i:') or tag.startswith('iq:i:') or tag.startswith('sq:i:'):
                        tag = tag[:A_end+5]
                        new_tags.append(tag)
                    else:
                        new_tags.append(tag)
                d['tags'] = new_tags
                x = pysam.AlignedSegment.from_dict(d, r.peer.header)
                f2.write(x)
        elif UMI_type == 'G5':
            G_start, G_end = find_Gstart(d['seq'])
            if G_start > 0:
                seq2 = d['seq'][:G_start]  # should be just UMI

                diff = len(seq2) - umi_len
                if diff < 0:  # UMI may have ended with Gs
                    seq2 = d['seq'][:G_start-diff]

                seq_extra = 'NA'
                if diff > 0:
                    seq_extra = seq2[:diff]
                    seq2 = seq2[diff:]

                rec = {'id': r.peer.qname,
                       'clip_len': len(seq2),
                       'extra': seq_extra,
                       'UMI': seq2,
                       'BC': 'NA',  # Brendan's current design has only UMI, no BC
                       'BC_rev': 'NA',
                       'BC_match': 'NA',
                       'BC_top_rank': 'NA'}
                writer1.writerow(rec)

                # subset the sequence to remove the UMIs and "G"s
                d['seq'] = d['seq'][G_end:]
                d['qual'] = d['qual'][G_end:]
                assert len(d['seq'])==len(d['qual'])
                new_tags = []
                for tag in d['tags']:
                    if tag.startswith('zs:B'):  # defunct CCS tag, don't use
                        pass
                    elif tag.startswith('dq:i:') or tag.startswith('iq:i:') or tag.startswith('sq:i:'):
                        tag = tag[:5] + tag[5+G_end:]
                        new_tags.append(tag)
                    else:
                        new_tags.append(tag)
                d['tags'] = new_tags
                x = pysam.AlignedSegment.from_dict(d, r.peer.header)
                f2.write(x)

    f1.close()
    f2.close()


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("bam_filename", help="CCS BAM with cDNA primer removed (post LIMA)")
    parser.add_argument("output_prefix", help="Output prefix")
    parser.add_argument("-u", "--umi_len", type=int, help="Length of UMI")
    parser.add_argument("-b", "--bc_len", type=int, help="Length of cell barcode")
    parser.add_argument("--umi_type", choices=['A3', 'G5'], help="Location of the UMI")
    parser.add_argument("--bc_rank_file", help="(Optional) cell barcode rank file from short read data")


    args = parser.parse_args()

    # ToDo: figure out later how to do top ranked barcodes for 10X data
    shortread_bc = {}  # dict of cell barcode -> "Y" for top ranked
    if args.bc_rank_file is not None:
        reader = DictReader(open(args.bc_rank_file), delimiter='\t')
        for r in reader:
            shortread_bc[r['cell_barcode']] = r['top_ranked']

    clip_out(args.bam_filename, args.umi_len, args.bc_len, args.output_prefix, args.umi_type, shortread_bc)
