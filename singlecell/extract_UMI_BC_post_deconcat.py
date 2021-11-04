#!/usr/bin/env python
import os, sys
from Bio.Seq import Seq
from csv import DictReader, DictWriter
import pysam
import parasail
import pdb

VALID_CIGAR_SYMBOL = ['I', 'D']

SCOREMAT = parasail.matrix_create("ACGT", 2, -5)


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

def find_G5_seq(query, G5_sequence, min_score):
    """
    Search for the 10X linker sequence (ex: TTTCTTATATGGG) and return the beginning, end of it
    ex:
    5' primer - (16bp BC) - (10bp UMI) - G5_sequence
    :return: (0-based start, 1-based end) or (-1, -1) if not found
    """
    o = parasail.sg_qx_trace(query, G5_sequence, 3, 1, SCOREMAT)
    if o.score >= min_score:
        return o.get_traceback().comp.find('|'), o.end_query+1
    else:
        return -1, -1


def clip_out(bam_filename, umi_len, bc_len, output_prefix, UMI_type, shortread_bc={},
             tso_len=0, g5_clip_seq=None, g3_clip_seq=None):
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
    assert UMI_type in ('G5-10X', 'A3-10X')
    umi_bc_len = umi_len + bc_len

    if UMI_type == 'G5-clip':
        try:
            import parasail
        except ImportError:
            print("need parasail library for G5-clip mode! Abort!", file=sys.stderr)
            sys.exit(-1)
        para_mat = parasail.matrix_create("ACGT", 2, -5)
        para_search_len = umi_len + len(g5_clip_seq) + 10

    FIELDS = ['id', 'old_id', 'clip_len', 'extra', 'UMI', 'BC', 'BC_rev', 'BC_match', 'BC_top_rank']
    if tso_len > 0: FIELDS += ['TSO']

    f1 = open(output_prefix + '.tagged.csv', 'w')
    writer1 = DictWriter(f1, FIELDS, delimiter='\t')
    writer1.writeheader()

    reader = pysam.AlignmentFile(open(bam_filename), 'rb', check_sq=False)
    #reader = BamIO.IndexedBamReader(bam_filename)
    #reader2 = pysam.AlignmentFile(open('flt.bam'), 'rb', check_sq=False)
    # mock the RG
    MOCK_RG_NAME = '38988035/1--0'
    header_dict = reader.header.as_dict()
    header_dict['RG'][0]['ID'] = MOCK_RG_NAME
    f2 = pysam.AlignmentFile(output_prefix+'.tagged.bam', 'wb', header=pysam.AlignmentHeader.from_dict(header_dict))

    counter = 0
    for r in reader:
        counter += 1
        if counter % 100000 == 0:
            print("Processed {0} reads...".format(counter))
        #if counter >= 50000: break
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


        if UMI_type == 'A3-10X':
            # (polyA tail) - (10bp UMI) - (16bp BC) - (possible leftover 3bp padding of AGA)
            assert g3_clip_seq is not None
            pad_len = len(g3_clip_seq)
            #pdb.set_trace()
            A_start, A_end = find_Aend(d['seq'])
            if A_end > 0:
                seq2 = d['seq'][A_end:]  # should be just UMI + BC, unless UMI started with 'A's

                diff = len(seq2) - umi_bc_len
                if diff < 0: # UMI may have started with 'A's, no padding
                    seq2 = d['seq'][A_end+diff:]
                elif diff == 0: # no padding, just UMI and BC
                    pass
                else: # some padding after UMI and BC. find and trim padding
                    found_pad_flag = False
                    pad_offset = pad_len
                    while pad_offset > 0:
                        if d['seq'][-pad_offset:] == g3_clip_seq[:pad_offset]:
                            found_pad_flag = True
                            break
                        pad_offset -= 1
                    if found_pad_flag:
                        seq2 = d['seq'][A_end:-pad_offset]
                        diff = len(seq2) - umi_bc_len
                        seq2 = d['seq'][A_end+diff:]

                seq_umi = seq2[:umi_len]
                seq_bc = seq2[umi_len:umi_bc_len]

                # reverse complement BC because it's always listed in rev comp in short read data
                seq_bc_rev = str(Seq(seq_bc).reverse_complement())

                match = 'Y' if seq_bc_rev in shortread_bc else 'N'
                match_top = 'Y' if (match=='Y' and shortread_bc[seq_bc_rev]=='Y') else 'N'

                d['name'] = d['name'].split('/')[0] + '/' + str(counter) + '/ccs'
                rec = {'id': d['name'],
                       'old_id': r.qname,
                       'clip_len': len(seq2),
                       'extra': 'NA',
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
                    if tag.startswith('zs:B'):  # defunct CCS tag, don't use
                        pass
                    elif tag.startswith('RG:Z:'): # ignore read group, put in a hard-coded one is isoseq3 refine is happy
                        new_tags.append('RG:Z:' + MOCK_RG_NAME)
                    elif tag.startswith('ql:') or tag.startswith('qt:') or tag.startswith('bl:') or tag.startswith('bt:'):
                        pass
                    elif tag.startswith('zm:i'):
                        new_tags.append('zm:i:' + str(counter)) # maunally change zmw becuz duplicated
                    elif tag.startswith('dq:i:') or tag.startswith('iq:i:') or tag.startswith('sq:i:'):
                        tag = tag[:A_end+5]
                        new_tags.append(tag)
                    else:
                        new_tags.append(tag)
                # adding UMI (XM) and BC (XC)
                new_tags.append('XM:Z:' + str(seq_umi))
                new_tags.append('XC:Z:' + str(seq_bc))
                new_tags.append('XA:Z:XM-XC')
                new_tags.append('XG:Z:GGG')
                d['tags'] = new_tags
                #d['name'] = d['name'].split('/')[0] + '/' + str(counter) + '/ccs'
                if len(d['seq']) >= 100 and len(seq_umi)==umi_len and len(seq_bc)==bc_len:
                    r_header_dict = r.header.as_dict()
                    r_header_dict['RG'][0]['ID'] = MOCK_RG_NAME
                    x = pysam.AlignedSegment.from_dict(d, pysam.AlignmentHeader.from_dict(r_header_dict))
                    f2.write(x)
                else:
                    print("Ignore {0} because: seqlen {1}, umilen {2}, bclen {3}".format(r.qname, len(d['seq']), len(seq_umi), len(seq_bc)))

        elif UMI_type == 'G5-10X':
            assert g5_clip_seq is not None
            assert tso_len == len(g5_clip_seq)

            umi_bc_tso_len = bc_len + umi_len + tso_len
            tso_min_score = (tso_len-2)*2
            G_start, G_end = find_G5_seq(d['seq'][:umi_bc_tso_len*2], g5_clip_seq, min_score=tso_min_score)

            if G_start >= 0:
                seq2 = d['seq'][:G_end]  # this is BC - UMI - TSO
                seq_tso = d['seq'][G_start:G_end]

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
                    #pdb.set_trace()
                    seq_extra = 'NA'
                    seq_bc = seq2[:bc_len+diff]
                    seq_umi = seq2[bc_len+diff:umi_bc_len+diff]

                # reverse complement BC because it's always listed in rev comp in short read data
                seq_bc_rev = str(Seq(seq_bc).reverse_complement())
                match = 'Y' if seq_bc_rev in shortread_bc else 'N'
                match_top = 'Y' if (match=='Y' and shortread_bc[seq_bc_rev]=='Y') else 'N'

                d['name'] = d['name'].split('/')[0] + '/' + str(counter) + '/ccs'
                rec = {'id': d['name'],
                       'old_id': r.qname,
                       'clip_len': len(seq2),
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
                    elif tag.startswith('RG:Z:'): # ignore read group, put in a hard-coded one is isoseq3 refine is happy
                        new_tags.append('RG:Z:' + MOCK_RG_NAME)
                    elif tag.startswith('ql:') or tag.startswith('qt:') or tag.startswith('bl:') or tag.startswith('bt:'):
                        pass
                    elif tag.startswith('zm:i'):
                        new_tags.append('zm:i:' + str(counter)) # maunally change zmw becuz duplicated
                    elif tag.startswith('dq:i:') or tag.startswith('iq:i:') or tag.startswith('sq:i:'):
                        tag = tag[:5] + tag[5+G_end:]
                        new_tags.append(tag)
                    else:
                        new_tags.append(tag)
                # adding UMI (XM) and BC (XC)
                new_tags.append('XM:Z:' + str(seq_umi))
                new_tags.append('XC:Z:' + str(seq_bc))
                new_tags.append('XA:Z:XM-XC')
                new_tags.append('XG:Z:GGG')
                d['tags'] = new_tags
                #d['name'] = d['name'].split('/')[0] + '/' + str(counter) + '/ccs'
                if len(d['seq']) >= 100 and len(seq_umi)==umi_len and len(seq_bc)==bc_len:
                    r_header_dict = r.header.as_dict()
                    r_header_dict['RG'][0]['ID'] = MOCK_RG_NAME
                    x = pysam.AlignedSegment.from_dict(d, pysam.AlignmentHeader.from_dict(r_header_dict))
                    f2.write(x)
                else:
                    print("Ignore {0} because: seqlen {1}, umilen {2}, bclen {3}".format(r.qname, len(d['seq']), len(seq_umi), len(seq_bc)))

    f1.close()
    f2.close()


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("bam_filename", help="CCS BAM after deconcation")
    parser.add_argument("output_prefix", help="Output prefix")
    parser.add_argument("-u", "--umi_len", type=int, help="Length of UMI")
    parser.add_argument("-b", "--bc_len", type=int, help="Length of cell barcode")
    parser.add_argument("-t", "--tso_len", type=int, default=0, help="Length of TSO (for G5-10X only)")
    parser.add_argument("--umi_type", choices=['A3-10X', 'G5-10X'], help="Location of the UMI")
    parser.add_argument('--g5_clip_seq', help="Sequence before UMI for G5-clip (for G5-10X only)")
    parser.add_argument("--a3_clip_seq", help="Extra padding sequence at end of A3 end (for A3-10X only)")
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

    clip_out(args.bam_filename, args.umi_len, args.bc_len, args.output_prefix, args.umi_type,
             shortread_bc, args.tso_len, args.g5_clip_seq, args.a3_clip_seq)
