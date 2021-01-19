#!/usr/bin/env python
import os, sys
from Bio.Seq import Seq
from csv import DictReader, DictWriter
import pysam
import pdb

VALID_CIGAR_SYMBOL = ['I', 'D']

def iter_cigar_string(cigar_string):
    # ex: 1=1X2=1X1=1D5=8D3=27I
    num = cigar_string[0]
    for s in cigar_string[1:]:
        if str.isdigit(s):
            num += s
        else:
            yield int(num), s
            num = ''

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
    if i >= 0:
        # now find the end
        j = i + 3
        while j < len_seq and seq[j]=='G':
            j += 1
        return i, j
    else:
        return -1, -1



def clip_out(bam_filename, umi_len, bc_len, output_prefix, UMI_type, shortread_bc={}, tso_len=0, g5_clip_seq=None):
    """
    :param bam_filename: BAM of post-LIMA (primer-trimmed) CCS sequences
    :param UMI_type: either 'A3' or 'G5' or 'G5-10X'
    :param shortread_bc: a dict of barcode -> "Y|N" for top-ranked. If given, came from short read data.

    --------
    G5-10X
    --------
    5' primer -- BC --- UMI -- TSO --- GGG --- transcript --- polyA

    --------
    G5-clip
    assumes input is like below, where the 5'/3' primer already removed by lima
    Here, we will only clip out the UMI, and write out the rest of the sequence, keeping the RT + transcript
    There is no assumption about the polyA tail existing or not
    --------
    5' primer -- UMI -- [RT primer] --- transcript --- 3' primer
    """
    assert UMI_type in ('A3', 'G5', 'G5-10X', 'G5-clip')
    umi_bc_len = umi_len + bc_len

    if UMI_type == 'G5-clip':
        try:
            import parasail
        except ImportError:
            print("need parasail library for G5-clip mode! Abort!", file=sys.stderr)
            sys.exit(-1)
        para_mat = parasail.matrix_create("ACGT", 2, -5)
        para_search_len = umi_len + len(g5_clip_seq) + 10

    FIELDS = ['id', 'clip_len', 'extra', 'UMI', 'BC', 'BC_rev', 'BC_match', 'BC_top_rank']
    if tso_len > 0: FIELDS += ['TSO']

    f1 = open(output_prefix + '.trimmed.csv', 'w')
    writer1 = DictWriter(f1, FIELDS, delimiter='\t', dialect='unix')
    writer1.writeheader()

    reader = pysam.AlignmentFile(bam_filename, 'rb', check_sq=False)
    #reader = BamIO.IndexedBamReader(bam_filename)
    f2 = pysam.AlignmentFile(output_prefix+'.trimmed.bam', 'wb', header=reader.header)

    for r in reader:
        d = r.to_dict()

        #is_rev_strand = r.flag >> 4 & 1
        if (r.flag >> 4 & 1):
            d['seq'] = str(Seq(r.seq).reverse_complement())
            d['qual'] = r.qual[::-1]
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

                if bc_len == 0:
                    seq_bc = ''
                else:
                    seq_bc = seq2[-bc_len:]

                if umi_len == 0:
                    seq_umi = ''
                else:
                    if bc_len == 0:
                        seq_umi = seq2[-umi_len:]
                    else:
                        seq_umi = seq2[-(bc_len+umi_len):-bc_len]


                # reverse complement BC because it's always listed in rev comp in short read data
                seq_bc_rev = str(Seq(seq_bc).reverse_complement())

                match = 'Y' if seq_bc_rev in shortread_bc else 'N'
                match_top = 'Y' if (match=='Y' and shortread_bc[seq_bc_rev]=='Y') else 'N'

                rec = {'id': r.qname,
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
                x = pysam.AlignedSegment.from_dict(d, r.header)
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

                rec = {'id': r.qname,
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
                x = pysam.AlignedSegment.from_dict(d, r.header)
                f2.write(x)
        elif UMI_type == 'G5-clip':
            o1 = parasail.sg_qx_trace(d['seq'][:para_search_len],
                                               g5_clip_seq,
                                               10,
                                               3,
                                               para_mat)

            #  'tags': ['bx:B:i,22,20',
            #   ...
            #   'qe:i:2835',
            #   'bc:B:S,0,1',
            #   'bl:Z:CCCGCGTGGCCTCCTGAATTAT',
            #   'bt:Z:CATTGCCACTGTCTTCTGCT',
            #   'RG:Z:70de1488/0--1']}
            c_num, c_type = next(iter_cigar_string(str(o1.cigar.decode, 'utf-8')))
            if c_type == 'I': # this is the (extra) + UMI
                seq2 = d['seq'][:c_num]
                seq_extra = 'NA'
                diff = len(seq2) - umi_len
                if diff < 0: # we need to get a few more bases from the primers
                    tag_dict = dict(x.split(':', 1) for x in d['tags'])
                    try:
                        if tag_dict['bc'] == 'B:S,0,1': # + strand
                            assert tag_dict['bl'].startswith('Z:')
                            Fseq = tag_dict['bl'][2:]  # trimming away the Z:
                        elif tag_dict['bc'] == 'B:S,1,0': # - strand
                            assert tag_dict['bt'].startswith('Z:')
                            Fseq = str(Seq(tag_dict['bt'][2:]).reverse_complement())
                        seq2 = Fseq[diff:] + seq2 # rescue bases from the trimmed F primer
                    except KeyError:
                        pass # just silently not do anything and output the shorter UMI
                        #print("WARNING: older version of lima output, lacking 'bc' tag. Ignoring read {0}...".format(r.qname))
                elif diff > 0:  # there's extras
                    seq_extra = seq2[:diff]
                    seq2 = seq2[diff:]

                rec = {'id': r.qname,
                       'clip_len': len(seq2),
                       'extra': seq_extra,
                       'UMI': seq2,
                       'BC': 'NA',  # Brendan's current design has only UMI, no BC
                       'BC_rev': 'NA',
                       'BC_match': 'NA',
                       'BC_top_rank': 'NA'}
                writer1.writerow(rec)

                # subset the sequence to remove the UMI (but keep the G5 clip seq)
                d['seq'] = d['seq'][c_num:]
                d['qual'] = d['qual'][c_num:]
                assert len(d['seq'])==len(d['qual'])
                new_tags = []
                for tag in d['tags']:
                    if tag.startswith('zs:B'):  # defunct CCS tag, don't use
                        pass
                    elif tag.startswith('dq:i:') or tag.startswith('iq:i:') or tag.startswith('sq:i:'):
                        tag = tag[:5] + tag[5+c_num:]
                        new_tags.append(tag)
                    else:
                        new_tags.append(tag)
                d['tags'] = new_tags
                x = pysam.AlignedSegment.from_dict(d, r.header)
                f2.write(x)
        elif UMI_type == 'G5-10X':
            # need to first invert the sequence so polyA is at the end
            d['seq'] = str(Seq(d['seq']).reverse_complement())
            d['qual'] = d['qual'][::-1]
            # now it is BC -- UMI -- TSO -- GGG -- transcript -- polyA
            umi_bc_tso_len = bc_len + umi_len + tso_len
            G_start, G_end = find_Gstart(d['seq'][umi_bc_tso_len:umi_bc_tso_len+10])

            #pdb.set_trace()

            if G_start >= 0:
                G_start += umi_bc_tso_len
                G_end   += umi_bc_tso_len

                seq2 = d['seq'][:G_start]  # this is BC - UMI - TSO
                seq_tso = seq2[-tso_len:] + d['seq'][G_start:G_end]

                diff = len(seq2) - umi_bc_tso_len
                if diff > 0: # beginning may have included untrimmed primers
                    seq_extra = seq2[:diff]
                    seq2 = seq2[diff:]
                    seq_bc = seq2[:bc_len]
                    seq_umi = seq2[bc_len:umi_bc_len]
                elif diff == 0:
                    seq_extra = 'NA'
                    seq_bc = seq2[:bc_len]
                    seq_umi = seq2[bc_len:umi_bc_len]
                elif diff < 0: # we may have accidentally trimmed away some bases for BC, can't do anything
                    seq_extra = 'NA'
                    seq_bc = seq2[:bc_len+diff]
                    seq_umi = seq2[bc_len+diff:umi_bc_len+diff]

                # reverse complement BC because it's always listed in rev comp in short read data
                seq_bc_rev = str(Seq(seq_bc).reverse_complement())
                match = 'Y' if seq_bc_rev in shortread_bc else 'N'
                match_top = 'Y' if (match=='Y' and shortread_bc[seq_bc_rev]=='Y') else 'N'

                rec = {'id': r.qname,
                       'clip_len': len(seq2)+(G_end-G_start),
                       'extra': seq_extra,
                       'UMI': seq_umi,
                       'BC': seq_bc,
                       'TSO': seq_tso,
                       'BC_rev': seq_bc_rev,
                       'BC_match': match,
                       'BC_top_rank': match_top}
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
                x = pysam.AlignedSegment.from_dict(d, r.header)
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
    parser.add_argument("-t", "--tso_len", type=int, default=0, help="Length of TSO (for G5-10X only)")
    parser.add_argument("--umi_type", choices=['A3', 'G5', 'G5-10X', 'G5-clip'], help="Location of the UMI")
    parser.add_argument('--g5_clip_seq', help="Sequence before UMI for G5-clip (for G5-clip only)")
    parser.add_argument("--bc_rank_file", help="(Optional) cell barcode rank file from short read data")


    args = parser.parse_args()

    if args.bc_len < 0:
        print("bc_len can't be a negative number!", file=sys.stderr)
        sys.exit(-1)
    if args.umi_len < 0:
        print("umi_len can't be a negative number!", file=sys.stderr)
        sys.exit(-1)
    if args.umi_len + args.bc_len <= 0:
        print("umi_len + bc_len must be at least 1 bp long!", file=sys.stderr)
        sys.exit(-1)

    # ToDo: figure out later how to do top ranked barcodes for 10X data
    shortread_bc = {}  # dict of cell barcode -> "Y" for top ranked
    if args.bc_rank_file is not None:
        reader = DictReader(open(args.bc_rank_file), delimiter='\t')
        for r in reader:
            shortread_bc[r['cell_barcode']] = r['top_ranked']

    clip_out(args.bam_filename, args.umi_len, args.bc_len, args.output_prefix, args.umi_type, shortread_bc, args.tso_len, args.g5_clip_seq)
