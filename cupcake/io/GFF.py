import os, re, sys
import pdb
from bx.intervals.intersection import Interval
from bx.intervals.intersection import IntervalNode
from bx.intervals.intersection import IntervalTree
from collections import defaultdict
from csv import DictReader

class GTF:
    def __init__(self, gtf_filename):
        self.gtf_filename = gtf_filename
        self.genome = defaultdict(lambda: IntervalTree()) # chr --> IntervalTree --> (0-start, 1-end, transcript ID)
        self.transcript = defaultdict(lambda: IntervalTree()) # tID --> IntervalTree --> (0-start, 1-end, {'ith': i-th exon, 'eID': exon ID})
        self.exon = defaultdict(lambda: []) # (0start,1end) --> list of (tID, ith-exon, chr)
        self.transcript_info = {} # tID --> chr
        
        self.readGTF(self.gtf_filename)
    
    def readGTF(self, filename):
        """
        GTF files
        (0) chr
        (1) annotation source
        (2) type: gene|transcript|CDS|exon|UTR
        (3) 1-based start 
        (4) 1-based end
        (5) ignore
        (6) strand: +|-
        (7) phase
        (8) extra stuff (gene ID, transcript ID...) 
        """
        for line in open(filename):
            if line.startswith('#'): continue # header section, ignore
            if len(line.strip()) == 0: continue # some gtf files have blank lines
            raw = line.strip().split('\t')
            chr = raw[0]
            type = raw[2]
            strand = raw[6]
            start0, end1 = int(raw[3])-1, int(raw[4])
            gtype, gstat = 'NA', 'NA'
            gName = 'NA'
            tName = 'NA'
            tSupportLevel = 'NA'
            gtags = []
            for stuff in raw[8].split('; '):
                _a, _b = stuff.split(None, 1)
                if _a == "transcript_id": tID = _b[1:-1] # removing quotes ""
                elif _a == "transcript_name": tName = _b[1:-1] # removing quotes
                elif _a == "transcript_support_level": tSupportLevel = _b[1:-1]
                elif _a == 'gene_name': gName = _b[1:-1]
                elif _a == "gene_id": gID = _b[1:-1] # removing quotes ""
                elif _a == "gene_type": gtype = _b[1:-1]
                elif _a == "gene_status": gstat = _b[1:-1]           
                elif _a == 'tag': gtags.append(_b[1:-1]) 
                
            if type == 'transcript':
                self.genome[chr].insert(start0, end1, tID)
                self.transcript_info[tID] = {'chr':chr, 'gname':gName, 'gid':gID, 'type': gtype, 'status': gstat, 'strand': strand, 'tags': gtags, 'tname':tName, 'tlevel': tSupportLevel}
                ith = 0
            elif type == 'exon':
                self.transcript[tID].insert(start0, end1, {'ith':ith,'chr':chr})
                self.exon[(start0,end1)].append((tID, ith, chr))
                ith += 1

        
    def get_exons(self, tID):
        """
        Return the list of intervals for a given tID
        """
        pp = []
        self.transcript[tID].traverse(pp.append)
        return pp
    
    def find(self, chr, start0, end1):
        return list(set(self.genome[chr].find(start0, end1)))

class polyAGFF(GTF):
    def readGTF(self, filename):
        with open(filename) as f:
            for line in f:
                if line.startswith('##'): continue # just comments
                raw = line.strip().split('\t')
                assert raw[2] in ("polyA_signal", "polyA_site", "pseudo_polyA")
                if raw[2] == "polyA_site":
                    chrom = raw[0]
                    strand = raw[6]
                    start = int(raw[3])-1
                    end = int(raw[4])
                    for stuff in raw[8].split('; '):
                        _a, _b = stuff.split(None, 1)
                        if _a == "transcript_id": tID = _b[1:-1] # removing quotes ""
                    self.genome[chrom].insert(start, end, tID)
                    self.transcript_info[tID] = {'chr':chrom, 'strand':strand}
                    self.transcript[tID].insert(start, end, 0)

class TSSGFF(GTF):
    def readGTF(self, filename):
        with open(filename) as f:
            for line in f:
                if line.startswith('##'): continue
                raw = line.strip().split('\t')
                assert raw[2] in ('Gencode TSS')
                chrom = raw[0]
                strand = raw[6]
                start = int(raw[3])-1
                end = int(raw[4])
                for stuff in raw[8].split('; '):
                    _a, _b = stuff.split(None, 1)
                    if _a == "gene_id": gID = _b[1:-1] # removing quotes ""
                self.genome[chrom].insert(start, end, (gID,start))
                assert start+1 == end
                if gID in self.transcript_info:
                    self.transcript_info[gID].append({'chr':chrom, 'strand':strand, 'position':start})
                else:
                    self.transcript_info[gID] = [{'chr':chrom, 'strand':strand, 'position':start}]
    
class ucscGTF(GTF):
    """
    UCSC-style GFF, which is
    0) seqname (chromosome)
    1) source
    2) feature (gene|exon|mRNA...)
    3) start (1-based)
    4) end (1-based)
    5) score
    6) strand
    7) frame 
    8) group
    """
    def readGTF(self, filename):
        for line in open(filename):
            raw = line.strip().split('\t')
            if raw[2] == 'exon':
                _chr = raw[0]
                strand = raw[6]
                start = int(raw[3])-1
                end = int(raw[4])
                tID = raw[8]
                
                if tID not in self.transcript:
                    # new transcript
                    self.genome[_chr].insert(start, end, tID)
                    self.transcript_info[tID] = {'chr':_chr, 'gid':None, 'type': None, 'status':None, 'strand': strand}
                    ith = 0
                    self.exon[(start,end)].append((tID, ith, _chr))
                    self.transcript[tID].insert(start, end, {'ith':ith,'chr':_chr})
                else:
                    ith += 1
                    self.transcript[tID].insert(start, end, {'ith':ith,'chr':_chr})
                    self.exon[(start,end)].append((tID, ith, _chr))
                    

class variantRecord:
    def __init__(self, chrom, type, start, end, reference, variant, freq, coverage, confidence):
        self.chr = chrom
        self.type = type
        self.start = start
        self.end = end
        self.reference = reference
        self.variant = variant
        self.freq = freq
        self.coverage = coverage
        self.confidence = confidence

    def __str__(self):
        return """
        {t} {c}:{s}-{e}
        reference: {ref}
        variant: {var} ({freq})
        coverage: {cov}
        confidence: {cof}
        """.format(t=self.type, c=self.chr, s=self.start, e=self.end, \
                   ref=self.reference, var=self.variant, freq=self.freq, \
                   cov=self.coverage, cof=self.confidence)


class variantGFFReader:
    """
    Chr1    .       substitution    86591   86591   .       .       .       reference=T;variantSeq=T/G;frequency=35/15;
coverage=66;confidence=40
    """
    def __init__(self, filename):
        self.filename = filename
        self.f = open(filename)
        while True:
            cur = self.f.tell()
            if not self.f.readline().startswith('#'): break
        self.f.seek(cur)

    def __iter__(self):
        return self

    def next(self):
        return self.read()

    def read(self):
        cur = self.f.tell()
        line = self.f.readline().strip()
        if self.f.tell() == cur:
            raise StopIteration("EOF reached!!")
        raw = line.strip().split('\t')

        chrom = raw[0]
        type = raw[2]
        start = int(raw[3])
        end = int(raw[4])
        for x in raw[8].split(';'):
            a, b = x.split('=')
            if a == 'reference': reference = b
            elif a == 'variantSeq': variant = b
            elif a == 'frequency': freq = b
            elif a == 'coverage': coverage = int(b)
            elif a == 'confidence': confidence = int(b)

        return variantRecord(chrom, type, start, end, reference, variant, freq, coverage, confidence)




class Coords(GTF):
    def readGTF(self, filename):
        """
        .coords files
        (0) gene name
        (1) chr
        (2) number of exons
        (3) strand
        (4) list of space-separated 1-based start, 1-based end
        """
        for line in open(filename):
            raw = line.strip().split()
            tID = raw[0]
            chr = raw[1]
            ith = 0
            
            if tID in self.transcript:
                print("duplicate tID {0} seen, ignore!".format(tID), file=sys.stderr)
                continue
            
            self.transcript_info[tID] = {'chr':chr}
            
            for i in xrange(4, len(raw), 2):
                start0 = int(raw[i])-1
                end1 = int(raw[i+1])
                self.genome[chr].insert(start0, end1, tID)
                self.transcript[tID].insert(start0, end1, {'ith':ith, 'chr':chr})
                self.exon[(start0, end1)].append((tID, ith, chr))
                
                i += 1
                
def write_gtf_records(gtf, tIDs, output_filename):
    f = open(output_filename, 'w')
    for tID in tIDs:
        info = gtf.transcript_info[tID]
        _chr = info['chr']
        strand = info['strand']
        f.write("{chr}\tJUNK\tgene\t{start}\t{end}\t.\t{strand}\t.\tgene_id \"{gid}\"; transcript_id \"{tid}\"")
        
    f.close()
    
        
            
class btabReader:
    def __init__(self, filename):
        self.filename = filename
        self.f = open(filename)

    def __iter__(self):
        return self
    
    def next(self):
        return self.read()            
            
    def read(self):
        """
        (0) chr
        (1)-(2) blank
        (3) gmap
        (4) blank
        (5) seqid
        (6) ref start (1-based)
        (7) ref end
        (8) seq start (1-based)
        (9) seq end
        (10) score
        (11)-(12) blank
        (13) ith-seq
        (14) ith-exon
        
        start-end will be flipped if start > end !!!
        """
        cur = self.f.tell()
        line = self.f.readline().strip()
        if cur == self.f.tell():
            raise StopIteration("EOF reached!")
        raw = line.split('\t')
        chr= raw[0]
        seqid = raw[5]
        rStart1 = int(raw[6])
        rEnd1 = int(raw[7])
        i = raw[-2]
        if rStart1 > rEnd1: rStart1, rEnd1 = rEnd1, rStart1
        return {'chr': chr, 'seqid': seqid, 'rStart1': rStart1, 'rEnd1': rEnd1, 'i': i}

class btabBlockReader(btabReader):
    def next(self):
        recs = [self.read()]
        while recs[-1]['i']==recs[0]['i']:
            cur = self.f.tell()
            recs.append(self.read())
        self.f.seek(cur)    
        return recs[:-1]

pbid_rex = re.compile('(PB\.\d+|PBfusion\.\d+)(\.\d+){0,1}')

class gmapRecord:
    def __init__(self, chr, coverage, identity, strand, seqid, geneid=None):
        """
        Record keeping for GMAP output:
        chr, coverage, identity, seqid, exons
        
        exons --- list of Interval, 0-based start, 1-based end
        """
        assert strand == '+' or strand == '-'
        self.chr = chr
        self.coverage = coverage
        self.identity = identity
        self.strand = strand
        self.seqid = seqid
        self.ref_exons = []
        self.seq_exons = []
        self.cds_exons = []
        self.scores = []

        #pdb.set_trace()
        self.geneid = None
        # handle gene ids specially for PB.X.Y and PBfusion.X.Y
        if geneid is not None:
            self.geneid = geneid
        else:
            m = pbid_rex.match(seqid)
            if m is not None:
                self.geneid = m.group(1)

    def __eq__(self, other):
        return self.chr==other.chr and \
               self.coverage == other.coverage and \
               self.identity == other.identity and \
               self.strand == other.strand and \
               self.seqid == other.seqid and \
               self.ref_exons == other.ref_exons and \
               self.seq_exons == other.seq_exons and \
               self.cds_exons  == other.cds_exons and \
               self.scores == other.scores and \
               self.geneid == other.geneid

    def __str__(self):
        return """
        chr: {0}
        strand: {1}
        coverage: {2}
        identity: {3}
        seqid: {4}
        geneid: {8}
        ref exons: {5}
        seq exons: {6}
        scores: {7}
        """.format(self.chr, self.strand, self.coverage, self.identity, self.seqid, self.ref_exons, self.seq_exons, self.scores, self.geneid)
        
    def __getattr__(self, key):
        if key == 'rstart' or key == 'start':
            return self.get_start()
        elif key == 'rend' or key == 'end':
            return self.get_end()
        else:
            raise AttributeError(key)
        
    def get_start(self): return self.ref_exons[0].start
    
    def get_end(self): return self.ref_exons[-1].end
        

    def add_cds_exon(self, start, end):
        self.cds_exons.append(Interval(start, end))


    def add_exon(self, rStart0, rEnd1, sStart0, sEnd1, rstrand, score):
        assert rStart0 < rEnd1 and sStart0 < sEnd1
        if rstrand == '-':
            assert len(self.ref_exons) == 0 or self.ref_exons[0].start >= rEnd1
            self.scores.insert(0, score)
            self.ref_exons.insert(0, Interval(rStart0, rEnd1))
        else:
            assert len(self.ref_exons) == 0 or self.ref_exons[-1].end <= rStart0
            self.scores.append(score)
            self.ref_exons.append(Interval(rStart0, rEnd1))
        if rstrand == '-':
            self.seq_exons.insert(0, Interval(sStart0, sEnd1))
        else:
            self.seq_exons.append(Interval(sStart0, sEnd1))
            
    
class gmapGFFReader(object):
    def __init__(self, filename):
        self.filename = filename
        self.f = open(filename)
        # read through possible header
        while True:
            cur = self.f.tell()
            if not self.f.readline().startswith('#') or self.f.tell()==cur: # first non-# seen or EOF
                self.f.seek(cur)
                break
        self.sanity_format_check()
        
    def __iter__(self):
        return self
    
    def __next__(self):
        return self.read()

    def sanity_format_check(self):
        """
        GFF3 formats are supposed to be tab-delimited and 9 required fields
        https://learn.gencore.bio.nyu.edu/ngs-file-formats/gff3-format/
        """
        cur = self.f.tell()
        raw = self.f.readline().strip().split('\t')
        if len(raw) != 9:
            print("ERROR:Sanity checking {0} is GFF3 format: expected 9 tab-delimited fields but saw {1}! Abort!".format(\
                self.filename, len(raw)))
            sys.exit(-1)
        self.f.seek(cur)

            
    def read(self):
        """
        GFF files
        (0) chr
        (1) annotation source
        (2) type: gene|transcript|CDS|exon|UTR
        (3) 1-based start # MUST CONVERT TO 0-based!!!
        (4) 1-based end
        (5) score (I think it's similarity for GMAP)
        (6) strand: +|-
        (7) phase
        (8) extra stuff (gene ID, transcript ID...) 
        
        For gmap output, a series is delimited by '###' line
        """
        cur = self.f.tell()
        line = self.f.readline().strip()
        if self.f.tell() == cur:
            raise StopIteration("EOF reached!!")
        raw = line.strip().split('\t')
        while raw[0].startswith('#'):
            line = self.f.readline().strip()
            raw = line.strip().split('\t')
        
        if len(raw) == 0 or raw[0]=='': 
            raise StopIteration("EOF reached!!")

        assert raw[2] == 'gene'
        raw = self.f.readline().strip().split('\t')
        assert raw[2] == 'mRNA'
        chr = raw[0]
        strand = raw[6]
        for blob in raw[8].split(';') :
            if blob.startswith('coverage='): coverage = float(blob[9:])
            elif blob.startswith('identity='): identity = float(blob[9:])
            elif blob.startswith('Name='): seqid = blob[5:]
      
        rec = gmapRecord(chr, coverage, identity, strand, seqid)
        
        cds_exons = []
        cds_seq_start = None
        cds_seq_end = None
        while True:
            line = self.f.readline().strip()
            if line.startswith('##'):
                rec.cds_exons = cds_exons
                rec.cds_seq_start = cds_seq_start
                rec.cds_seq_end = cds_seq_end
                return rec
            raw = line.split('\t')
            type = raw[2]            
            if type == 'exon':
                rstart1, rend1 = int(raw[3]), int(raw[4])
                score = float(raw[5])
                rstrand = raw[6] # this is the strand on the reference genome
                for blob in raw[8].split(';'):
                    if blob.startswith('Target='):
                        # sstrand is the strand on the query sequence
                        junk, sstart1, send1, sstrand = blob.split()
                        sstart1 = int(sstart1)
                        send1 = int(send1)
                        rec.sstrand = sstrand
                try:
                    rec.add_exon(rstart1-1, rend1, sstart1-1, send1, rstrand, score)
                except AssertionError:
                    print("{0} has non-colinear exons!".format(rec.seqid), file=sys.stderr)
                    while True:
                        line = self.f.readline().strip()
                        if line.startswith('##'): return rec
                rec.strand = rstrand
            elif type == 'CDS':
                rstart1, rend1 = int(raw[3]), int(raw[4])
                cds_exons.append(Interval(rstart1-1, rend1))
                for blob in raw[8].split(';'):
                    if blob.startswith('Target='):
                        junk, sstart1, send1, sstrand = blob.split()
                        sstart1 = int(sstart1)
                        send1 = int(send1)
                        cds_seq_start = sstart1-1 if cds_seq_start is None else cds_seq_start
                        cds_seq_end = send1
            else:
                raise Exception("Not supposed to see type {0} here!!".format(type))

        return rec
    
            
class pasaGFFReader(gmapGFFReader):
    """
    Slight differences in PASA's GTF output (.gtf)
    Each transcript is separated by 1 or more blank lines
    """
    def read(self):
        cur = self.f.tell()
        line = self.f.readline().strip()
        if self.f.tell() == cur:
            raise StopIteration("EOF reached!!")
        while line.startswith('#'): # header section, ignore
            line = self.f.readline().strip()              
        raw = line.split('\t')
        assert raw[2] == 'transcript'
        
        chr = raw[0]
        strand = raw[6]
        for blob in raw[8].split('; '):
            if blob.startswith('transcript_id'): # ex: transcript_id "asmbl_7"
                tid = blob[15:-1] 
            elif blob.startswith('gene_id'): # ex: gene_id "S2"
                gid = blob[9:-1]
    
        rec = gmapRecord(chr=chr, coverage=None, identity=None, strand=strand, seqid=tid, geneid=gid)
        
        while True:
            #pdb.set_trace()
            line = self.f.readline().strip()
            if line.startswith('###'): # end of this record
                return rec                
            raw = line.split('\t')
            type = raw[2]
            start1, end1 = int(raw[3]), int(raw[4])
            if type == 'exon':
                rec.add_exon(start1-1, end1, -2, -1, None)

def write_collapseGFF_format(f, r):
    f.write("{chr}\tPacBio\ttranscript\t{s}\t{e}\t.\t{strand}\t.\ttranscript_id \"{tid}\"; gene_id \"{gid}\";\n".format(chr=r.chr, s=r.start+1, e=r.end, strand=r.strand,gid=r.geneid, tid=r.seqid))
    for exon in r.ref_exons:
        f.write("{chr}\tPacBio\texon\t{s}\t{e}\t.\t{strand}\t.\ttranscript_id \"{tid}\"; gene_id \"{gid}\";\n".format(chr=r.chr, s=exon.start+1, e=exon.end, strand=r.strand, gid=r.geneid, tid=r.seqid))
    if r.cds_exons is not None:
        for exon in r.cds_exons:
            f.write("{chr}\tPacBio\tCDS\t{s}\t{e}\t.\t{strand}\t.\ttranscript_id \"{tid}\"; gene_id \"{gid}\";\n".format(chr=r.chr, s=exon.start+1, e=exon.end, strand=r.strand, gid=r.geneid, tid=r.seqid))


class collapseGFFReader(gmapGFFReader):
    def read(self):
        """
        PacBio-style GFF from the collapsed output, which is
        0) chrmosome
        1) source (PacBio)
        2) feature (transcript|exon)
        3) start (1-based)
        4) end (1-based)
        5) score (always .)
        6) strand
        7) frame (always .)
        8) blurb

        ex: 
        chr1    PacBio  transcript      897326  901092  .       +       .       gene_id "PB.1"; transcript_id "PB.1.1";
        chr1    PacBio  exon    897326  897427  .       +       .       gene_id "PB.1"; transcript_id "PB.1.1";
        """
        cur = self.f.tell()
        line = self.f.readline().strip()
        if self.f.tell() == cur:
            raise StopIteration("EOF reached!!")
                
        raw = line.strip().split('\t')
        assert raw[2] == 'transcript'
        chr = raw[0]
        strand = raw[6] 
        seqid = None
        geneid = None
        for stuff in raw[8].split(';'):
            if len(stuff.strip()) > 0:
                a, b = stuff.strip().split()
                if a == 'transcript_id':
                    seqid = b[1:-1]
                if a == 'gene_id':
                    geneid = b[1:-1]

        rec = gmapRecord(chr, coverage=None, identity=None, strand=strand, seqid=seqid, geneid=geneid)
        
        while True:
            cur = self.f.tell()
            line = self.f.readline().strip()
            if self.f.tell() == cur:
                return rec
            raw = line.split('\t')
            if raw[2] == 'exon':
                s, e = int(raw[3])-1, int(raw[4])
                rec.add_exon(s, e, s, e, rstrand='+', score=None)
            elif raw[2] == 'CDS':
                s, e = int(raw[3])-1, int(raw[4])
                rec.add_cds_exon(s, e)
            else: # another new record, wind back and return
                self.f.seek(cur)
                return rec
        raise Exception("Should not reach here!")


fusion_seqid_rex = re.compile('(\S+\\.\d+)\\.(\d+)')
class collapseGFFFusionReader(collapseGFFReader):
    def read(self):
        r0 = super(collapseGFFFusionReader, self).read()
        m = fusion_seqid_rex.match(r0.seqid)
        fusion_id = m.group(1)
        records = [r0]
        while True:
            cur = self.f.tell()
            try:
                r = super(collapseGFFFusionReader, self).read()
            except StopIteration:
                return fusion_id, records
            m = fusion_seqid_rex.match(r.seqid)
            if m.group(1) != fusion_id:
                self.f.seek(cur)
                return fusion_id, records
            else:
                records.append(r)




class ucscGFFReader(gmapGFFReader):
    def read(self):
        """
        UCSC-style GFF, which is
        0) seqname (chromosome)
        1) source
        2) feature (gene|exon|mRNA...)
        3) start (1-based)
        4) end (1-based)
        5) score
        6) strand
        7) frame 
        8) group
        
        A series is delimited by '###' line
        """
        cur = self.f.tell()
        line = self.f.readline().strip()
        if self.f.tell() == cur:
            raise StopIteration("EOF reached!!")
                
        raw = line.strip().split('\t')
        assert raw[2] == 'exon'
        chr = raw[0]
        s, e = int(raw[3])-1, int(raw[4])
        strand = raw[6] 
        seqid = raw[8]       
      
        rec = gmapRecord(chr, coverage=None, identity=None, strand=strand, seqid=seqid)
        rec.add_exon(s, e, s, e, strand, score=None)
        
        while True:
            line = self.f.readline().strip()
            if line.startswith('###'):
                return rec
            raw = line.split('\t')
            assert raw[2] == 'exon'
            s, e = int(raw[3])-1, int(raw[4]) 
            rec.add_exon(s, e, s, e, strand, score=None)
        return rec        
                
def GFFReader(filename):
    """
    Reads the 2nd column to decide GMAP or PASA parser to use
    """
    with open(filename) as f:
        program = f.readline().strip().split('\t')[1]
    if program == 'GMAP':
        return gmapGFFReader(filename)
    elif program == 'PASA':
        return pasaGFFReader(filename)
    else:
        raise Exception("{0} is not a recognizable GFF program".format(program))

def write_fancyGeneformat(f, r):
    for exon in r.ref_exons:
        f.write("{0} exon {1} {2}\n".format(r.seqid, exon.start+1, exon.end+1))

def write_GFF_UCSCformat(f, r):  
    """
    UCSC GTF format:
    0) seqname
    1) source
    2) feature (gene|exon|mRNA...)
    3) start (1-based)
    4) end (1-based)
    5) score
    6) strand
    7) frame 
    8) group
    
    r should be gmapRecord object
    """  
    ref_exons = r.ref_exons
    if r.strand == '-':
        ref_exons.reverse()
    for exon in r.ref_exons:
        f.write(r.chr + '\t')
        f.write('NA\t')
        f.write("exon\t")
        f.write(str(exon.start+1) +'\t')
        f.write(str(exon.end) + '\t')
        f.write(".\t")
        f.write(r.strand + '\t')
        try:
            f.write(str(r.score)+'\t')
        except:
            f.write(".\t")
        f.write(r.seqid + '\n')
    f.write('###\n')  
    
def convert_BLAST9rec_to_gmapRecord(rec_list):
    """
    Adds .chr, .seqid, and .ref_exons so we can use it to write in UCSC format
    """
    assert len(rec_list) > 0
    chr = rec_list[0].sID
    seqid = rec_list[0].qID
    strand = rec_list[0].strand
    assert all(x.sID==chr for x in rec_list)
    assert all(x.qID==seqid for x in rec_list)
    assert all(x.strand==strand for x in rec_list)
    
    r = gmapRecord(chr, coverage=0, identity=0, strand=strand, seqid=seqid)
    r.ref_exons = [Interval(x.sStart, x.sEnd) for x in rec_list]
    
    return r
              
def btab_reclist_to_interval_list_0basedStart(recs):
    """
    Return chr, list of IntervalNode
    """
    tree = IntervalTree()
    for rec in recs:
        tree.insert(rec['rStart1']-1, rec['rEnd1'])
    path = []
    tree.traverse(path.append)
    chr = recs[0]['chr']
    return chr, path

def getOverlap(a, b):
    return max(0, min(a.end, b.end) - max(a.start, b.start))
    
  
def CompareSimCoordinatesToAlnPath(alnPath, simCoordinates):
    #
    # do silly little dynamic programming to align sets of exons.
    # This could be done in a while loop if there is a 1-1
    # correspondende of exons that overlap, but if multiple overlap,
    # that could cause problems.
    # 
    nAlnExons = len(alnPath)
    nSimExons = len(simCoordinates)
    scoreMat = [[0 for j in xrange(nSimExons+1) ] for i in xrange(nAlnExons+1) ]
    pathMat  = [[0 for j in xrange(nSimExons+1) ] for i in xrange(nAlnExons+1) ]

    diagonal = 0
    up = 1
    left = 2

    for i in xrange(nAlnExons):
        pathMat[i+1][0] = up
    for j in xrange(nSimExons):
        pathMat[0][j+1] = left
    pathMat[0][0] = diagonal
    #return 0

    for i in xrange(nAlnExons):
        for j in xrange(nSimExons):
            overlapScore = 0
            if len(simCoordinates[j].find(alnPath[i].start, alnPath[i].end)) > 0: # overlaps!
                overlapScore = getOverlap(alnPath[i], simCoordinates[j])*1./(simCoordinates[j].end-simCoordinates[j].start) # GetOverlapPercent(alnPair, simCoordinates.exonList[j])
                assert 0 <= overlapScore <= 1.
                scoreMat[i+1][j+1] = scoreMat[i][j] + overlapScore
                pathMat[i+1][j+1]  = diagonal
            else :
                order = simCoordinates[j].end <= alnPath[i].start #WhichIntervalIsFirst(alnPair, simCoordinates.exonList[j])
                if order:
                    scoreMat[i+1][j+1] = scoreMat[i][j+1] -2  # penalize assembled exons that were skipped
                    pathMat[i+1][j+1]  = up
                else:             
                    scoreMat[i+1][j+1] = scoreMat[i+1][j] -1 # penalize gencode exons being skipped
                    pathMat[i+1][j+1]  = left

    #pdb.set_trace()
    i = nAlnExons    
    j = nSimExons
    matchedExons = []
    _cur_best_j = nSimExons
    for j in xrange(nSimExons-1, -1, -1):
        if scoreMat[i][j] > scoreMat[i][_cur_best_j]:
            _cur_best_j = j
    j = _cur_best_j
    while (i > 0 and j > 0):
        if (pathMat[i][j] == diagonal):
            matchedExons.append((j-1,i-1)) # format should be (ref_ind, seq_ind)
            i = i - 1
            j = j - 1
        elif(pathMat[i][j] == left):
            j = j - 1
        else:
            i = i - 1
    matchedExons.reverse()
    return (scoreMat[nAlnExons][_cur_best_j]-(nSimExons-_cur_best_j), matchedExons)  
        
        
def match_transcript(gtf, chr, exon_path):
    """
    exon_tree is an IntervalTree, so it's already sorted
    """
    num_exon = len(exon_path)
    
    #print 'matching transcript for:', exon_path
    
    best_score, best_matchedExons, best_tID, best_tNum = 0, None, None, None
    for tID in gtf.find(chr, exon_path[0].start, exon_path[-1].end):
        t_paths = gtf.get_exons(tID)
        
        score, matchedExons = CompareSimCoordinatesToAlnPath(exon_path, t_paths)
       
        #print 'matching:', tID, score, matchedExons
        #pdb.set_trace()
        if score > best_score:
            best_tID = tID
            best_tNum = len(t_paths)
            best_score = score
            best_matchedExons = matchedExons
            
    return {'score':best_score, 'matchedExons':best_matchedExons, 'tID': best_tID, 'tID_num_exons': best_tNum}


def categorize_transcript_recovery(info):
    """
    full --- means that every exon in the tID was covered!
    fused --- full, but assembled exon match start > 0, meaning
              likely fusion of overlapped transcripts
    5missX --- means that the assembled one is missing beginning X exons
    3missY --- means that the assembled one is missing ending Y exons
    skipped --- means that the asseembled one is missing some intermediate exons!
    """
    if len(info['matchedExons']) == info['tID_num_exons']: 
        if info['matchedExons'][0][1] == 0: return 'full'
        else: return 'fused'
    msg = ''
    if info['matchedExons'][0][0] > 0: 
        msg += '5miss' if info['strand'] == '+' else '3miss'
        msg += str(info['matchedExons'][0][0])
    if info['matchedExons'][-1][0] < info['tID_num_exons']-1: 
        msg += (';' if msg!='' else '') 
        msg += '3miss' if info['strand'] == '+' else '5miss' 
        msg += str(info['tID_num_exons']-1-info['matchedExons'][-1][0])
        
    if msg == '': # must be missing some ground truth exons!
        return 'skipped'
    return msg

def evaluate_alignment_boundary_goodness(ref_exons, aln_exons, matches):
    """
    Returns a list of comma-separated numbers (head,tail).
    For each head element: 0 if precise, +k if seq starts at ref.start+k, -k if ref.start-k
    For each tail element: 0 if precise, +k if seq starts at ref.end+k, -k if ref.end-k
    """
    result = []
    for ind_ref, ind_aln in matches:
        result.append((aln_exons[ind_aln].start-ref_exons[ind_ref].start, \
                       aln_exons[ind_aln].end-ref_exons[ind_ref].end))
    return result
        
        
    
    

def main(gtf):
    transcript_tally = {}
    for tID in gtf.transcript: 
        transcript_tally[tID] = [0]*len(gtf.get_exons(tID))
    for r in btabBlockReader('sim_gencode_20x_first1000_test2.gmap.tophits.btab'):
        path = btab_reclist_to_interval_list(r)
        info = match_transcript(gtf, r[0]['chr'], path)
        if info['matchedExons'] is None:
            print("Did not find a match for {0}!".format(r[0]['seqid']), file=sys.stderr)
            continue
        for i, j in info['matchedExons']:
            transcript_tally[info['tID']][i] += 1
    return transcript_tally
    
def main_pasa(gtf):
    pasa_tally = {}
    for tID in gtf.transcript:
        pasa_tally[tID] = [0]*len(gtf.get_exons(tID))
    pasa = GTF('sim_gencode_20x_first1000_test2.pasa_assemblies.denovo_transcript_isoforms.gtf')
    for tID in pasa.transcript:
        path = pasa.get_exons(tID)
        chr = pasa.exon[(path[0].start,path[0].end)][0][2]
        
        info = match_transcript(gtf, chr, path)
        if info['matchedExons'] is None:
            print("Did not find a match for {0}!".format(tID), file=sys.stderr)
            continue
        for i, j in info['matchedExons']:
            pasa_tally[info['tID']][i] += 1
    return pasa_tally


def eval_gmap(gtf, gmap_filename, input_filename):
    """
    Expected seqID format: m000000_000000_00000_cSIMULATED_s0_p0/0/0_1250 or p0/1395/ccs
    
    Input: 
    gtf --- GTF/Coords object as ground truth transcripts
    gmap_filename --- gmap output in .gff format
    input_filename --- input fasta to gmap (to identify unmapped seqs)
    
    Output: <output_prefix> is just <gmap_filename>
    <output_prefix>.bad --- list of seqids that had no GMAP output or did not match a transcript
    <output_prefix>.report --- 
      <seqid>, <seqlen>, <seqMatchStart>, <seqMatchEnd>, <transcript/gene ID>, <category:full|5missX|3missY|skipped>, <matchedExons>
    """
    from Bio import SeqIO
    output_prefix = gmap_filename
    fbad = open(output_prefix+'.bad', 'w')
    fbad.write("seqID\tinfo\n")
    fgood = open(output_prefix+'.report', 'w')
    fgood.write("seqID\tseqLen\tchr\tstrand\tseqMatchStart0\tseqMatchEnd1\trefID\tcategory\tmatches\tboundary\n")
    
    seqlen_dict = dict((r.id, len(r.seq)) for r in SeqIO.parse(open(input_filename),'fastq' if input_filename.endswith('.fastq') else 'fasta'))
    seqid_missed = seqlen_dict.keys()
    
    for rec in gmapGFFReader(gmap_filename):
        chr = rec.chr        
        seqid = rec.seqid
        print("seqid: {0}".format(seqid))
        seqlen = seqlen_dict[seqid]
        try:
            seqid_missed.remove(seqid)
        except ValueError: # already removed, ignore?
            pass
        info = match_transcript(gtf, chr, rec.ref_exons)
        info['strand'] = rec.strand
        if info['matchedExons'] is None:
            fbad.write("{0}\tBAD\n".format(seqid))
        else:
            fgood.write("{seqid}\t{seqlen}\t{chr}\t{strand}\t{smstart0}\t{smend1}\t{refID}\t{cat}\t{mat}\t{bound}\n".format(\
                seqid=seqid, seqlen=seqlen, chr=chr, smstart0=rec.start, smend1=rec.end,\
                strand=rec.strand,\
                refID=info['tID'], cat=categorize_transcript_recovery(info), mat=info['matchedExons'],\
                bound=evaluate_alignment_boundary_goodness(gtf.get_exons(info['tID']), rec.ref_exons, info['matchedExons'])))
        
    for seqid in seqid_missed:
        fbad.write("{0}\tMISSED\n".format(seqid))
    fbad.close()
    fgood.close()
    

def eval_pasa(gtf, pasa_filename, gmap_report_filename):
    """
    
    Output:
    <gID> <tID> <number of exons> <refID> <category> <matches> 

    """
    output_prefix = pasa_filename
    fbad = open(output_prefix+'.bad', 'w')
    fbad.write("tID\tinfo\n")
    fgood = open(output_prefix+'.report', 'w')
    fgood.write("gID\ttID\tchr\tstrand\tnum_exon\ttlen\trefID\treflen\trefStrand\tcategory\tmatches\tboundary\n")
    
    refid_missed = list(set(x['refID'] for x in DictReader(open(gmap_report_filename), delimiter='\t')))

    for rec in pasaGFFReader(pasa_filename):
        gid = rec.seqid
        tid = rec.seqid
        num_exon = len(rec.ref_exons)
        tLen = sum(x.end-x.start for x in rec.ref_exons) # i know it's confusing but seq_exon is not used in parsing pASA output!

        info = match_transcript(gtf, rec.chr, rec.ref_exons)
        info['strand'] = rec.strand
        refid = info['tID']
        try:
            refid_missed.remove(refid)
        except ValueError:
            pass
        refLen = sum(x.end-x.start for x in gtf.get_exons(info['tID']))
        if info['matchedExons'] is None:
            fbad.write("{0}\tBAD\n".format(tid))
        else:
            fgood.write("{gid}\t{tid}\t".format(gid=gid, tid=tid))
            fgood.write("{chr}\t{strand}\t".format(chr=rec.chr, strand=rec.strand))
            fgood.write("{num_exon}\t{tLen}\t".format(num_exon=num_exon, tLen=tLen))
            fgood.write("{refID}\t{refLen}\t".format(refID=info['tID'], refLen=refLen))
            fgood.write("{refStrand}\t".format(refStrand=gtf.transcript_info[refid]['strand']))
            fgood.write("{cat}\t{mat}\t".format(cat=categorize_transcript_recovery(info), mat=info['matchedExons']))
            fgood.write("{bound}\n".format(bound=evaluate_alignment_boundary_goodness(gtf.get_exons(info['tID']), rec.ref_exons, info['matchedExons'])))

    for refid in refid_missed:
        fbad.write("{0}\tMISSED\n".format(refid))
    fbad.close()
    fgood.close()
    
                                                                      
def make_exon_report(gtf, gmap_report_filename):
    """
    Output for each exon:
    <tID>   <exon number 0-based>  <length>  <coverage>
    
    Output will be written to .exon_report   
    """     
    coverage = defaultdict(lambda: defaultdict(lambda: 0)) # tID --> ith-exon --> count           
    for r in DictReader(open(gmap_report_filename), delimiter='\t'):
        tID = r['refID']
        for i,j in eval(r['matches']):
            coverage[tID][i] += 1
    
    f = open(gmap_report_filename + '.exon_report', 'w')
    for tID in coverage:
        path = gtf.get_exons(tID)
        for ith, exon in enumerate(path):
            f.write("{0}\t{1}\t{2}\t{3}\n".format(tID, ith, exon.end-exon.start, coverage[tID][ith]))
                        
    f.close()
        
"""
##gff-version 3
# generated on Tue Dec 10 12:13:05 2013 by ./dump_all_gramene_dbs_continuing.pl
# for species maize
# genebuild 2010-01-MaizeSequence
5       ensembl gene    1579    3920    .       -       .       ID=GRMZM2G356204;Name=GRMZM2G356204;biotype=protein_coding
5       ensembl mRNA    1579    3920    .       -       .       ID=GRMZM2G356204_T01;Parent=GRMZM2G356204;Name=GRMZM2G356204_T01;biotype=protein_coding
5       ensembl exon    1579    3920    .       -       .       Parent=GRMZM2G356204_T01;Name=exon.12
5       ensembl CDS     1681    3903    .       -       .       Parent=GRMZM2G356204_T01;Name=CDS.13
5       ensembl gene    10731   23527   .       -       .       ID=GRMZM2G054378;Name=GRMZM2G054378;biotype=protein_coding
5       ensembl mRNA    22087   23527   .       -       .       ID=GRMZM2G054378_T09;Parent=GRMZM2G054378;Name=GRMZM2G054378_T09;biotype=protein_coding
5       ensembl intron  22956   23034   .       -       .       Parent=GRMZM2G054378_T09;Name=intron.19
5       ensembl intron  22799   22898   .       -       .       Parent=GRMZM2G054378_T09;Name=intron.20
5       ensembl intron  22634   22722   .       -       .       Parent=GRMZM2G054378_T09;Name=intron.21
5       ensembl intron  22456   22553   .       -       .       Parent=GRMZM2G054378_T09;Name=intron.22
5       ensembl exon    23035   23527   .       -       .       Parent=GRMZM2G054378_T09;Name=exon.23
5       ensembl exon    22899   22955   .       -       .       Parent=GRMZM2G054378_T09;Name=exon.24
5       ensembl exon    22723   22798   .       -       .       Parent=GRMZM2G054378_T09;Name=exon.25
5       ensembl exon    22554   22633   .       -       .       Parent=GRMZM2G054378_T09;Name=exon.26
5       ensembl exon    22087   22455   .       -       .       Parent=GRMZM2G054378_T09;Name=exon.27
5       ensembl CDS     23035   23193   .       -       .       Parent=GRMZM2G054378_T09;Name=CDS.28
5       ensembl CDS     22929   22955   .       -       0       Parent=GRMZM2G054378_T09;Name=CDS.29
"""
class MaizeGFFReader(collapseGFFReader):
    def read(self):
        """
        PacBio-style GFF from the collapsed output, which is
        0) chrmosome
        1) source (PacBio)
        2) feature (transcript|exon)
        3) start (1-based)
        4) end (1-based)
        5) score (always .)
        6) strand
        7) frame (always .)
        8) blurb

        ex:
        chr1    PacBio  transcript      897326  901092  .       +       .       gene_id "PB.1"; transcript_id "PB.1.1";
        chr1    PacBio  exon    897326  897427  .       +       .       gene_id "PB.1"; transcript_id "PB.1.1";
        """
        cur = self.f.tell()
        line = self.f.readline().strip()
        if self.f.tell() == cur:
            raise StopIteration("EOF reached!!")

        raw = line.strip().split('\t')
        if raw[2] == 'gene': # ignore this and read the next line 'mRNA'
            line = self.f.readline().strip()
            raw = line.strip().split('\t')
        assert raw[2] == 'mRNA'
        chr = raw[0]
        strand = raw[6]
        seqid = None
        for stuff in raw[8].split(';'): # ex: ID=GRMZM2G054378;Name=GRMZM2G054378;biotype=protein_coding
            a, b = stuff.strip().split('=')
            if a == 'ID':
                seqid = b
                break

        rec = gmapRecord(chr, coverage=None, identity=None, strand=strand, seqid=seqid)

        while True:
            cur = self.f.tell()
            line = self.f.readline().strip()
            if self.f.tell() == cur:
                return rec
            raw = line.split('\t')
            if raw[2] == 'exon':
                s, e = int(raw[3])-1, int(raw[4])
                rec.add_exon(s, e, s, e, rstrand=strand, score=None)
            elif raw[2] == 'CDS':
                s, e = int(raw[3])-1, int(raw[4])
                rec.add_cds_exon(s, e)
            elif raw[2] == 'intron':
                pass # ignore intron annotations
            else: # another new record, wind back and return
                self.f.seek(cur)
                return rec
        raise Exception("Should not reach here!")
                
    


"""
Because f***ing exonerate cannot generate GFF3 formats.

X       exonerate:est2genome    gene    78421768        78454596        5775    -       .       gene_id 0 ; sequence lcl|NM_001001171.1_cds_NP_001001171.1_1 ; gene_orient
ation + ; identity 99.92 ; similarity 99.92
X       exonerate:est2genome    utr5    78454482        78454596        .       -       .
X       exonerate:est2genome    exon    78454482        78454596        .       -       .       insertions 0 ; deletions 0 ; identity 100.00 ; similarity 100.00
X       exonerate:est2genome    splice5 78454480        78454481        .       -       .       intron_id 1 ; splice_site "GT"
X       exonerate:est2genome    intron  78451617        78454481        .       -       .       intron_id 1
X       exonerate:est2genome    splice3 78451617        78451618        .       -       .       intron_id 0 ; splice_site "AG"

"""
class ExonerateGFF2Reader(collapseGFFReader):
    def read(self):
        """
        PacBio-style GFF from the collapsed output, which is
        0) chrmosome
        1) source (PacBio)
        2) feature (transcript|exon)
        3) start (1-based)
        4) end (1-based)
        5) score (always .)
        6) strand
        7) frame (always .)
        8) blurb

        ex:
        chr1    PacBio  transcript      897326  901092  .       +       .       gene_id "PB.1"; transcript_id "PB.1.1";
        chr1    PacBio  exon    897326  897427  .       +       .       gene_id "PB.1"; transcript_id "PB.1.1";
        """
        cur = self.f.tell()
        line = self.f.readline().strip()
        if self.f.tell() == cur:
            raise StopIteration("EOF reached!!")

        raw = line.strip().split('\t')
        if raw[2] == 'gene': # read the gene line to get the ID and strand
            for stuff in raw[8].split(' ; '):
                a, b = stuff.strip().split(' ')
                if a == 'sequence':
                    seqid = b


        chr = raw[0]
        strand = raw[6]
        assert strand in ('+', '-')

        rec = gmapRecord(chr, coverage=None, identity=None, strand=strand, seqid=seqid)

        while True:
            cur = self.f.tell()
            line = self.f.readline().strip()
            if self.f.tell() == cur:
                return rec
            raw = line.split('\t')
            if raw[2] == 'exon':
                s, e = int(raw[3])-1, int(raw[4])
                rec.add_exon(s, e, s, e, rstrand=strand, score=None)
            elif raw[2] == 'CDS':
                pass # ignore CDS annotations
            elif raw[2] == 'intron':
                pass # ignore intron annotations
            elif raw[2] == 'similarity':
                # end of the record! return
                return rec
            elif raw[2] in ('splice3', 'splice5', 'utr3', 'utr5'):
                pass # ignore

        raise Exception("Should not reach here!")
            
        
        
        
        
        
    
    
