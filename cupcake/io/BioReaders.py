"""
Should always be faithful duplicate of sequence/BioReaders.py
Duplicated here for tofu installation. This one is called via cupcake.io.BioReaders.
"""

import re, sys, pdb
from collections import namedtuple
import pysam

Interval = namedtuple('Interval', ['start', 'end'])
                                 
class SimpleSAMReader:
    """
    A simplified SAM reader meant for speed. Skips CIGAR & FLAG parsing; identity/coverage calculation.
    """
    SAMheaders = ['@HD', '@SQ', '@RG', '@PG', '@CO']    
    def __init__(self, filename, has_header):
        self.filename = filename
        self.f = open(filename)
        self.header = ''
        if has_header:
            while True:
                cur = self.f.tell()
                line = self.f.readline()
                if line[:3] not in SimpleSAMReader.SAMheaders:
                    break
                self.header += line
            self.f.seek(cur)
    
    def __iter__(self):
        return self
    
    def __next__(self):
        line = self.f.readline().strip()
        if len(line) == 0:
            raise StopIteration
        return SimpleSAMRecord(line)    
    
  
class SimpleSAMRecord:
    cigar_rex = re.compile('(\d+)([MIDSHN])')
    SAMflag = namedtuple('SAMflag', ['is_paired', 'strand', 'PE_read_num', 'is_secondary', 'is_supplementary'])
    def __init__(self, record_line):
        """
        Simple bare bones version: only has
        
        qID, sID, sStart, sEnd, qStart, qEnd, cigar
        
        Simplified assumptions:
        -- must be end-to-end alignment (so qStart always 0)
        -- must be unspliced (no 'N' in cigar string)
        """
        self.qID = None
        self.sID = None
        self.sStart = None
        self.sEnd = None
        self.qStart = 0
        self.qEnd = None # length of SEQ
        self.cigar = None

        self.process(record_line)

    def __str__(self):
        msg = \
        """
        qID: {q}
        sID: {s}
        sStart-sEnd: {ss}-{se}
        qStart-qEnd: {qs}-{qe}
        cigar: {c}
        """.format(q=self.qID, s=self.sID, \
            ss=self.sStart, se=self.sEnd, qs=self.qStart, qe=self.qEnd, c=self.cigar)
        return msg

    def parse_cigar(self, cigar, start):
        """
        M - match
        I - insertion w.r.t. to ref
        D - deletion w.r.t. to ref
        N - skipped (which means splice junction)
        S - soft clipped
        H - hard clipped (not shown in SEQ)
        = - read match
        X - read mismatch

        ex: 50M43N3D

        NOTE: sets qStart & qEnd, which are often incorrect because of different ways to write CIGAR strings
              instead rely on XS/XE flags (from blasr or pbalign.py) to overwrite this later!!!

        Returns: genomic segment locations (using <start> as offset)
        """
        cur_end = start
        q_aln_len = 0
        for (num, type) in re.findall('(\d+)(\S)', cigar):
            num = int(num)
            if type == 'I':
                q_aln_len += num
            elif type in ('M', '=', 'X'):
                cur_end += num
                q_aln_len += num
            elif type == 'D':
                cur_end += num
        self.qEnd = self.qStart + q_aln_len
        self.sEnd = cur_end

            
    def process(self, record_line):
        """
        Only process cigar to get qEnd and sEnd
        """
        raw = record_line.split('\t')
        self.qID = raw[0]
        self.sID = raw[2]
        if self.sID == '*': # means no match! STOP here
            return
        self.sStart = int(raw[3]) - 1
        self.cigar = raw[5]
        self.parse_cigar(self.cigar, self.sStart)
        #self.flag = SimpleSAMRecord.parse_sam_flag(int(raw[1]))

    

class SAMReader:
    SAMheaders = ['@HD', '@SQ', '@RG', '@PG', '@CO']    
    def __init__(self, filename, has_header, ref_len_dict=None, query_len_dict=None):
        self.filename = filename
        self.f = open(filename)
        self.header = ''
        self.ref_len_dict = ref_len_dict
        self.query_len_dict = query_len_dict
        if has_header:
            while True:
                cur = self.f.tell()
                line = self.f.readline()
                if line[:3] not in SAMReader.SAMheaders:
                    break
                self.header += line
            self.f.seek(cur)
    
    def __iter__(self):
        return self
        
    def __next__(self):
        line = self.f.readline().strip()
        if len(line) == 0:
            raise StopIteration
        return SAMRecord(line, self.ref_len_dict, self.query_len_dict)        
    

class SAMRecord:
    SAMflag = namedtuple('SAMflag', ['is_paired', 'strand', 'PE_read_num', 'is_secondary', 'is_supplementary'])
    def __init__(self, record_line=None, ref_len_dict=None, query_len_dict=None):
        """
        Designed to handle BowTie SAM output for unaligned reads (PE read not yet supported)
        Can handle map to transfrag (no splicing) and genome (splicing)
        """
        self.qID = None
        self.sID = None
        self.sStart = None
        self.sEnd = None
        self.segments = None
        self.num_nonmatches = None
        self.num_ins = None
        self.num_del = None
        self.num_mat_or_sub = None

        self.qCoverage = None
        self.sCoverage = None

        self.sLen = None
        self.qLen = None
        self.cigar_qLen = None # qLen based on parsing CIGAR string
        # qStart, qEnd might get changed in parse_cigar
        self.qStart = 0
        self.qEnd = None # length of SEQ

        self.cigar = None
        self.flag = None

        self.identity = None
        self.record_line = record_line
        if record_line is not None:
            self.process(record_line, ref_len_dict, query_len_dict)

    def __str__(self):
        msg =\
        """
        qID: {q}
        sID: {s}
        cigar: {c}
        sStart-sEnd: {ss}-{se}
        qStart-qEnd: {qs}-{qe}
        segments: {seg}
        flag: {f}
        
        coverage (of query): {qcov}
        coverage (of subject): {scov}
        alignment identity: {iden}
        """.format(q=self.qID, s=self.sID, seg=self.segments, c=self.cigar, f=self.flag,\
            ss=self.sStart, se=self.sEnd, qs=self.qStart, qe=self.qEnd, iden=self.identity,\
            qcov=self.qCoverage, scov=self.sCoverage)
        return msg

    def __eq__(self, other):
        return self.qID == other.qID and self.sID == other.sID and\
               self.sStart == other.sStart and self.sEnd == other.sEnd and\
               self.segments == other.segments and self.qCoverage == other.qCoverage and\
               self.sCoverage == other.sCoverage and self.qLen == other.qLen and\
               self.sLen == other.sLen and self.qStart == other.qStart and\
               self.cigar == other.cigar and self.flag == other.flag and self.identity == other.identity

    @property
    def ref_exons(self):
        return self.segments

    def process(self, record_line, ref_len_dict, query_len_dict):
        """
        If SAM is from pbalign.py output, then have flags:
            XS: 1-based qStart, XE: 1-based qEnd, XQ: query length, NM: number of non-matches

        ignore_XQ should be False for BLASR/pbalign.py's SAM, True for GMAP's SAM
        
        0. qID
        1. flag
        2. sID
        3. 1-based offset sStart
        4. mapping quality (ignore)
        5. cigar
        6. name of ref of mate alignment (ignore)
        7. 1-based offset sStart of mate (ignore)
        8. inferred fragment length (ignore)
        9. sequence (ignore)
        10. read qual (ignore)
        11. optional fields
        """
        raw = record_line.split('\t')
        self.qID = raw[0]
        self.sID = raw[2]
        if self.sID == '*': # means no match! STOP here
            return
        self.sStart = int(raw[3]) - 1
        self.cigar = raw[5]
        self.segments = self.parse_cigar(self.cigar, self.sStart)
        self.sEnd = self.segments[-1].end
        self.flag = SAMRecord.parse_sam_flag(int(raw[1]))
        
        # process optional fields
        # XM: number of mismatches
        # NM: edit distance (sub/ins/del)
        for x in raw[11:]:
            if x.startswith('NM:i:'):
                self.num_nonmatches = int(x[5:])

        if ref_len_dict is not None:
            self.sCoverage = (self.sEnd - self.sStart) * 1. / ref_len_dict[self.sID]
            self.sLen = ref_len_dict[self.sID]

        if self.flag.strand == '-' and self.qLen is not None:
            self.qStart, self.qEnd = self.qLen - self.qEnd, self.qLen - self.qStart
         
        if query_len_dict is not None: # over write qLen and qCoverage, should be done LAST
            self.qLen = query_len_dict[self.qID]
            self.qCoverage = (self.qEnd - self.qStart) * 1. / self.qLen
            
        if self.num_nonmatches is not None:
            self.identity = 1. - (self.num_nonmatches * 1. / (self.num_del + self.num_ins + self.num_mat_or_sub))
            

    def parse_cigar(self, cigar, start):
        """
        M - match
        I - insertion w.r.t. to ref
        D - deletion w.r.t. to ref
        N - skipped (which means splice junction)
        S - soft clipped
        H - hard clipped (not shown in SEQ)
        = - read match
        X - read mismatch

        ex: 50M43N3D

        NOTE: sets qStart & qEnd, which are often incorrect because of different ways to write CIGAR strings

        Returns: genomic segment locations (using <start> as offset)
        """
        segments = []
        cur_start = start
        cur_end = start
        first_thing = True
        q_aln_len = 0
        self.num_del = 0
        self.num_ins = 0
        self.num_mat_or_sub = 0
        self.cigar_qLen = 0
        for (num, type) in re.findall('(\d+)(\S)', cigar):
            num = int(num)
            self.cigar_qLen += num
            if type == 'H' or type == 'S':
                if first_thing:
                    self.qStart += num
            elif type == 'I':
                q_aln_len += num
                self.num_ins += num
            elif type in ('M','=','X'):
                cur_end += num
                q_aln_len += num
                self.num_mat_or_sub += num
            elif type == 'D':
                cur_end += num
                self.num_del += num
            elif type == 'N': # junction, make a new segment
                segments.append(Interval(cur_start, cur_end))
                cur_start = cur_end + num
                cur_end = cur_start
            else:
                raise Exception("Unrecognized cigar character {0}!".format(type))
            first_thing = False
        if cur_start != cur_end:
            segments.append(Interval(cur_start, cur_end))
        self.qEnd = self.qStart + q_aln_len
        return segments

    @classmethod
    def parse_sam_flag(self, flag):
        """
		Heng Li's SAM https://samtools.github.io/hts-specs/SAMv1.pdf
        1 -- read is one of a pair
        2 -- alignment is one end of proper PE alignment          (IGNORE)
        4 -- read has no reported alignments                      (IGNORE)
        8 -- read is one of a pair and has no reported alignments (IGNORE)
        16 -- reverse ref strand
        32 -- other mate is aligned to ref strand
        64 -- first mate in pair
        128 -- second mate in pair
        256 -- not primary alignment
		512 -- not passing filters
		1024 -- PCR or optical duplicate
		2048 -- supplementary alignment

        Return: SAMflag
        """
        PE_read_num = 0
        strand = '+'
        is_supp = False
        is_secondary = False
        if flag >= 2048: # supplementary alignment
            flag -= 2048
            is_supp = True
        if flag >= 1024: #PCR or optical duplicate, should never see this...
            flag -= 1024
        if flag >= 512: #not passing QC, should never see this
            flag -= 512
        if flag >= 256: #secondary alignment, OK to see this if option given in BowTie
            flag -= 256
            is_secondary = True
        if flag >= 128:
            PE_read_num = 2
            flag -= 128
        elif flag >= 64:
            PE_read_num = 1
            flag -= 64
        if flag >= 32:
            flag -= 32
        if flag >= 16:
            strand = '-'
            flag -= 16
        if flag >= 8:
            flag -= 8
        if flag >= 4:
            flag -= 4
        if flag >= 2:
            flag -= 2
        assert flag == 0 or flag == 1
        is_paired = flag == 1
        return SAMRecord.SAMflag(is_paired, strand, PE_read_num, is_secondary, is_supp)
            

            
class GMAPSAMReader(SAMReader):
    def __next__(self):
        while True:
            line = self.f.readline().strip()
            if len(line) == 0:
                raise StopIteration
            if not line.startswith('@'): # header can occur at file end if the SAM was sorted
                break
        return GMAPSAMRecord(line, self.ref_len_dict, self.query_len_dict)
    
class GMAPSAMRecord(SAMRecord):
    def process(self, record_line, ref_len_dict=None, query_len_dict=None):
        """
        SAM files from pbalign.py have following optional fields:
            XS: 1-based qStart, XE: 1-based qEnd, XQ: query length, NM: number of non-matches
    
        0. qID
        1. flag
        2. sID
        3. 1-based offset sStart
        4. mapping quality (ignore)
        5. cigar
        6. name of ref of mate alignment (ignore)
        7. 1-based offset sStart of mate (ignore)
        8. inferred fragment length (ignore)
        9. sequence (ignore)
        10. read qual (ignore)
        11. optional fields
        """
        raw = record_line.split('\t')
        self.qID = raw[0]
        self.sID = raw[2]
        if self.sID == '*': # means no match! STOP here
            return
        self.sStart = int(raw[3]) - 1
        self.cigar = raw[5]
        self.segments = self.parse_cigar(self.cigar, self.sStart)
        self.sEnd = self.segments[-1].end
        self.flag = SAMRecord.parse_sam_flag(int(raw[1])) # strand can be overwritten by XS:A flag
        self._flag_strand = self.flag.strand # serve as backup for debugging

        for x in raw[11:]:
            if x.startswith('NM:i:'): # number of non-matches
                self.num_nonmatches = int(x[5:])
                self.identity = 1. - (self.num_nonmatches * 1. / (self.num_del + self.num_ins + self.num_mat_or_sub))
            elif x.startswith('XS:A:'): # strand information
                _s = x[5:]
                if _s!='?':
                    self._flag_strand = self.flag.strand # serve as backup for debugging
                    self.flag = SAMRecord.SAMflag(self.flag.is_paired, _s, self.flag.PE_read_num, is_secondary=False, is_supplementary=False)

        if ref_len_dict is not None:
            self.sCoverage = (self.sEnd - self.sStart) * 1. / ref_len_dict[self.sID]
            self.sLen = ref_len_dict[self.sID]

        if self.flag.strand == '-' and self.qLen is not None:
            self.qStart, self.qEnd = self.qLen - self.qEnd, self.qLen - self.qStart

        if self.qLen is not None:
            self.qCoverage = (self.qEnd - self.qStart) * 1. / self.qLen
           
        if query_len_dict is not None: # over write qLen and qCoverage, should be done LAST
            try:
                self.qLen = query_len_dict[self.qID]
            except KeyError: # HACK for blasr's extended qID
                raise Exception("Unable to find qID {0} in the input fasta/fastq!".format(self.qID))
            self.qCoverage = (self.qEnd - self.qStart) * 1. / self.qLen


class SplicedBAMReader:
    """
    The SplicedBAMReader imitates the behavior of GMAPSAMReader,
    basically accepted an aligned BAM file instead of aligned SAM file
    The returned records will have the same format as GMAPSAMRecord
    """

    def __init__(self, filename, ref_len_dict=None, query_len_dict=None):
        self.filename = filename
        self.reader = pysam.AlignmentFile(open(filename), 'rb', check_sq=False)
        self.ref_len_dict = ref_len_dict
        self.query_len_dict = query_len_dict

    def __iter__(self):
        return self

    def __next__(self):
        return self.grab_next_record()

    def grab_next_record(self):
        try:
            r = next(self.reader)
            samrec = SAMRecord(None) # we initiate the record and fill it in manually

            samrec.qID = r.qname
            samrec.qStart = r.qstart
            samrec.qEnd = r.qend
            if self.query_len_dict is not None:
                samrec.qLen = self.query_len_dict[samrec.qID]
            else:
                samrec.qLen = r.qlen
            samrec.sID = '*' if r.is_unmapped else r.reference_name
            samrec.sStart = r.reference_start
            samrec.sEnd = r.reference_end
            if self.ref_len_dict is not None:
                samrec.sLen = self.ref_len_dict[samrec.sID]
            else:
                samrec.sLen = r.reference_length
            samrec.cigar = r.cigarstring
            samrec.record_line = r.tostring()
            if samrec.sID == '*': # unmapped, nothing to parse
                return samrec

            # calling parse_cigar also sets num_ins, num_del, num_mat_or_sub, cigar_qlen
            samrec.segments = samrec.parse_cigar(r.cigarstring, r.reference_start)
            samrec.flag = SAMRecord.parse_sam_flag(r.flag)

            tag_d = dict(r.tags)
            if 'NM' in tag_d:   # this is used by regular minimap2
                samrec.num_nonmatches = tag_d['NM']
                samrec.identity = 1 - (samrec.num_nonmatches / samrec.qLen)
            else: # parse the cigar tuple to get the number of mismatches/insertions/deletions
                # https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.cigartuples
                # NOTE: we rely on this being run with mismatches represented by X
                samrec.num_nonmatches = 0
                for cigartype,cigarcount in r.cigartuples:
                    # 1:I, 2:D, 8:X
                    if cigartype in [1, 2, 8]:
                        samrec.num_nonmatches += cigarcount
                samrec.identity = 1 - (samrec.num_nonmatches / samrec.qLen)

            samrec.qCoverage = (r.qend-r.qstart)/samrec.qLen
            samrec.sCoverage = (r.reference_end-r.reference_start)/samrec.sLen
            return samrec
        except StopIteration:
            raise StopIteration


class SplicedBAMReaderRegioned(SplicedBAMReader):
    """
    Extension of SpliceBAMReader, except that an upfront [start_index, end_index) is defined,
    so that it'll return the records within that region

    The returned records will have the same format as GMAPSAMRecord
    """
    def __init__(self, filename, start_index, end_index, ref_len_dict=None, query_len_dict=None):
        self.filename = filename
        self.reader = pysam.AlignmentFile(open(filename), 'rb', check_sq=False)
        self.ref_len_dict = ref_len_dict
        self.query_len_dict = query_len_dict
        self.start_index = start_index
        self.end_index = end_index
        self.cur_index = start_index

        if self.start_index is None or self.end_index is None or self.start_index >= self.end_index:
            raise Exception(f"SplicedBAMReaderRegioned must be given proper integer [start_index, end_index)! Instead got {start_index}, {end_index}")

        for i in range(self.start_index):
            r = next(self.reader)

    def __iter__(self):
        return self

    def __next__(self):
        self.cur_index += 1
        if self.cur_index > self.end_index:
            raise StopIteration
        return self.grab_next_record()


