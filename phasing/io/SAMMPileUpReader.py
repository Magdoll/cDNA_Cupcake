__author__ = 'etseng@pacb.com'

"""
Parser for `samtools mpileup`

http://www.htslib.org/doc/samtools-1.1.html
1. chr
2. 1-based position
3. ref base
4. coverage
5. readBase
6. base qualities
7. alignment qualities

readBase:
.  match to ref
,  match to ref on rev
> or <    ref skipping  (ex: like 37N)
ACGTN  mismatch on + strand
acgn   mismatch on - strand
+{number}{AGCTNagctn} - insertion of some {number}
-{number}{...} deletion of some {number}  # also means in next {number}, you will see a *
^ begin of read, followed by asci-33 for quality
$ end of read
"""

import os, sys, re
import pdb
from collections import Counter

class MPileUpRecord(object):
    def __init__(self, chr, pos, ref, cov, readBase, baseQuals, alnQuals):
        """
        In addition to storing the 7 cols from mpileup,
        nalso stores
        counter: Counter of (key) -> (obs count in pileup)
        """
        self.chr = chr
        self.pos = pos
        self.ref = ref.upper() # let ref base always be upper case
        self.cov = cov
        self.nCov = None # this is the coverage of non-indel, non-skipped, which would be ACGTNacgtn
        self.nType = None # this is the number of non-indel, non-skipped bases accumulated at this record
        self.readBase = readBase
        self.baseQuals = baseQuals
        self.alnQuals = alnQuals

        self.counts = Counter()
        self.parse_readBase()

    def __str__(self):
        return """
        chr: {c}
        pos: {p} (1-based)
        ref: {r}
        cov: {v}
        nCov: {n}
        counts: {t}""".format(c=self.chr, p=self.pos+1, r=self.ref, v=self.cov, n=self.nCov, t=self.counts)

    def parse_readBase(self):
        """
        fill in self.counts
        """
        def not_indel_end_pos(i):
            return i >= len(self.readBase)-1 or self.readBase[i+1] not in ('+', '-', '$')

        rex = re.compile('(\d+)')
        def read_indel(start_index):
            m = rex.search(self.readBase, start_index)
            num = int(self.readBase[m.start():m.end()])
            return m.start(), m.end()+num

        sanity_counter = 0 # use this to track how many "reads" we've parsed to make sure parsing is correct
                           # this number should agree with self.cov which is 4-th column in mpileup
        i = 0 # pointer for current location in string self.readBase
        while i < len(self.readBase):
            b = self.readBase[i]
            if b in '<>': # ignore skipped refs
                sanity_counter += 1
                i += 1
                continue
            elif b == '*': # deletion, just advance
                i += 1
                sanity_counter += 1
                continue
            elif b == '^': # start of read followed by ascii and either a comma or dot (ex: ^I.)
                i += 3
                sanity_counter += 1
                continue
            elif b == '$': # end of read, DO NOT advance counter
                i += 1
                continue
            elif b == '.': # could be followed by indels or $, careful don't double count
                self.counts[self.ref] += 1
                sanity_counter += 1
                i += 1
            elif b == ',': # # could be followed by indels or $, careful don't double count
                self.counts[self.ref.lower()] += 1
                sanity_counter += 1
                i += 1
            elif b in 'ATCGNatcgn':
                self.counts[b] += 1
                sanity_counter += 1
                i += 1
            elif b == '-': # DO NOT ADVANCE the sanity counter! otherwise double counting
                start, end = read_indel(i+1)
                self.counts["-"+self.readBase[start:end]] += 1
                i = end
            elif b == '+': # insertion should be +{number}{bases}
                start, end = read_indel(i+1)
                self.counts["+"+self.readBase[start:end]] += 1
                i = end
            else:
                raise Exception("Unknown {0} in readBase!".format(b))

        assert self.cov == sanity_counter or (self.readBase=='*' and self.cov==0)
        # set nCov which is cov provided by non-indel non-skipped bases
        self.nCov = 0
        self.nType = 0
        for x in 'ATCGNatcgn':
            self.nCov += self.counts[x]
            if self.counts[x] > 0: self.nType += 1


class MPileUpReader(object):
    def __init__(self, filename):
        self.filename = filename
        self.f = open(filename)

    def __iter__(self):
        return self

    def __next__(self):
        cur = self.f.tell()
        line = self.f.readline()
        if self.f.tell() == cur:
            raise StopIteration
        return self.parseLine(line)

    def parseLine(self, line):
        raw = line.strip().split('\t')
        if (len(raw)==7 or len(raw)==15):
            cov = int(raw[3])
            #if cov > 0:
            return MPileUpRecord(chr=raw[0],\
                                pos=int(raw[1])-1,\
                                ref=raw[2],
                                cov=int(raw[3]),
                                readBase=raw[4],
                                baseQuals=raw[5],
                                alnQuals=raw[6])
        elif len(raw)==4:
            # only way to have only 4 columns is because after --min-BQ filtering there are no bases
            # ex:
            # fake    8728    T       3       .$.$.   ;q:     ]]]
            # fake    8729    T       0
            return MPileUpRecord(chr=raw[0],\
                                 pos=int(raw[1])-1,\
                                 ref=raw[2],
                                 cov=0,
                                 readBase='',
                                 baseQuals='',
                                 alnQuals='')
        else:
            raise Exception("Expected to have 7 cols in mpileup record \
            but saw only {0}, abort! Line was: {1}".format(len(raw), line))


