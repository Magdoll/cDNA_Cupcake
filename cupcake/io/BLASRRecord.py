__author__ = 'etseng@pacb.com'

#### THIS SHOULD BE A FAITHFUL REPLICATE OF BLASR IO from pbtranscript!


"""Define BLASR M4 reader, M5 reader."""


#M4DELIMITER, M5DELIMITER = " ", " "
# The M4 header including 13 fields.
# M4FIELDS = ["qName", "tName", "score", "percentSimilarity",
#            "qStrand", "qStart", "qEnd", "qLength",
#            "tStrand", "tStart", "tEnd", "tLength",
#            "mapQV"]
# M5FIELDS = ["qName", "qLength", "qStart", "qEnd", "qStrand",
#            "tName", "tLength", "tStart", "tEnd", "tStrand",
#            "score", "numMatch", "numMismatch", "numIns", "numDel", "mapQV",
#            "qAlignedSeq", "matchPattern", "tAlignedSeq"]
#
#M4HEADER = M4DELIMITER.join (M4FIELDS)
#M4FIELDSID = dict(zip(M4FIELDS, xrange(0, len(M4FIELDS))))
#


class BLASRRecord(object):

    """BLASR alignment."""

    def __init__(self, qID, qLength, qStart, qEnd, qStrand,
                 sID, sLength, sStart, sEnd, sStrand,
                 score, mapQV,
                 nMatch=None, nMismatch=None, nIns=None, nDel=None,
                 qAln=None, alnStr=None, sAln=None,
                 strand=None, identity=None):

        # define BLASR fields.
        self.qID = qID
        self.qLength = qLength
        self.qStart = qStart
        self.qEnd = qEnd

        def strandString(s):
            d = {'+': '+', '-': '-', 0: '+', 1: '-', '0': '+', '1': '-'}
            try:
                return d[s]
            except KeyError:
                raise ValueError("{s} is not a valid strand.".format(s=s))
        self.qStrand = strandString(qStrand)
        self.sID = sID
        self.sLength = sLength
        self.sStart = sStart
        self.sEnd = sEnd
        self.sStrand = strandString(sStrand)
        self.score = score
        self.nMatch = nMatch
        self.nMismatch = nMismatch
        self.nIns = nIns
        self.nDel = nDel
        self.mapQV = mapQV

        self.qAln = qAln
        self.alnStr = alnStr
        self.sAln = sAln

        self.strand = strandString(strand)
        self.identity = identity

    @property
    def alignment_length(self):
        """Alignment length."""
        return len(self.alnStr) if self.alnStr is not None \
            else None

    @property
    def mismatches(self):
        """Number of mismatches in an alignment."""
        return self.nMismatch

    @property
    def gaps(self):
        """Number of indels in an alignment."""
        return self.nIns + self.nDel \
            if self.nIns is not None and self.nDel is not None \
            else None

    @property
    def e(self):
        """NA"""
        return "NA"

    def __eq__(self, another):
        return self.__dict__ == another.__dict__

    def __str__(self):
        msg = """
        qID: {qID}
        sID: {sID}
        strand: {strand}
        identity: {identity}
        qMatch: {qStart}-{qEnd}
        sMatch: {sStart}-{sEnd}
        e-value: {e}
        BLASR score: {score}\n""".format(qID=self.qID, sID=self.sID,
                                         identity=self.identity,
                                         qStart=self.qStart, qEnd=self.qEnd,
                                         sStart=self.sStart, sEnd=self.sEnd,
                                         e=self.e, score=self.score, strand=self.strand)
        padding = " " * 8
        if self.alignment_length is not None:
            msg += "{p}alignment length: {alignment_length}\n".format(
                p=padding, alignment_length=self.alignment_length)
        if None not in [self.nMatch, self.nMismatch,
                        self.nIns, self.nDel, self.gaps]:
            msg += "{p}mismatches: {mismatches}\n{p}gaps: {gaps}\n".format(
                mismatches=self.nMismatch, gaps=self.gaps, p=padding)
        return msg

    @classmethod
    def fromM4String(cls, line, delimiter=' ', flipQS=False):
        """Interprets a string as BLASR M4 record."""
        try:
            fields = line.rstrip().split(delimiter)
            #d = dict(zip(fields, M4FIELDS))
            assert(len(fields) == 13)
            # trim qName of the last /0_len which is added by BLASR
            qName = fields[0]
            # if qName.rfind('/') > 0:
            #    qName = qName[:qName.rfind('/')]
            if qName.count('/') > 2:
                qName = "/".join(qName.split("/")[0:3])
            print qName
            tName = fields[1]
            score = int(fields[2])
            identity = float(fields[3])
            qStrand = fields[4]
            qStart = int(fields[5])
            qEnd = int(fields[6])
            qLength = int(fields[7])
            tStrand = fields[8]
            tStart = int(fields[9])
            tEnd = int(fields[10])
            tLength = int(fields[11])
            mapQV = int(fields[12])
            strand = '+' if qStrand == tStrand else '-'

            if flipQS:  # query is S, target is Q
                return BLASRRecord(
                    qID=tName, qLength=tLength,
                    qStart=tStart, qEnd=tEnd, qStrand=tStrand,
                    sID=qName, sLength=qLength,
                    sStart=qStart, sEnd=qEnd, sStrand=qStrand,
                    score=score, mapQV=mapQV,
                    strand=strand, identity=identity)
            else:  # query is Q, target is S
                return BLASRRecord(
                    qID=qName, qLength=qLength,
                    qStart=qStart, qEnd=qEnd, qStrand=qStrand,
                    sID=tName, sLength=tLength,
                    sStart=tStart, sEnd=tEnd, sStrand=tStrand,
                    score=score, mapQV=mapQV,
                    strand=strand, identity=identity)
        except (AssertionError, ValueError):
            errMsg = "String not recognized as a valid BLASR M4 record."
            raise ValueError(errMsg)

    @classmethod
    def fromM5String(cls, line):
        """Interprets a line as a BLASR M5 record."""
        try:
            fields = line.strip().split()
            assert(len(fields) == 19)
            qID = fields[0]
            # if qID.rfind('/') > 0:
            #    qID = qID[:qID.rfind('/')]
            if qID.count('/') > 2:
                qID = '/'.join(qID.split('/')[0:3])
                # blasr will add /<start>_<end> remove it
            qLength = int(fields[1])
            qStart = int(fields[2])
            qEnd = int(fields[3])
            qStrand = fields[4]
            sID = fields[5]
            sLength = int(fields[6])
            sStart = int(fields[7])
            sEnd = int(fields[8])
            sStrand = fields[9]
            score = int(fields[10])
            nMatch = int(fields[11])
            nMismatch = int(fields[12])
            nIns = int(fields[13])
            nDel = int(fields[14])
            mapQV = int(fields[15])
            qAln = fields[16]
            alnStr = fields[17]
            sAln = fields[18]
            #qCoverage = (qEnd - qStart) * 1. / qLength
            #sCoverage = (sEnd - sStart) * 1. / sLength
            # This identity is slightly different from percentSimilarity
            # of BLASR (percentSimilarity= (nMatch * 2 * 100 / len(q) + len(s))
            identity = float("{0:.4f}".format(nMatch * 100. / len(alnStr)))
            strand = '+' if qStrand == sStrand else '-'
            return BLASRRecord(
                qID=qID, qLength=qLength,
                qStart=qStart, qEnd=qEnd, qStrand=qStrand,
                sID=sID, sLength=sLength,
                sStart=sStart, sEnd=sEnd, sStrand=sStrand,
                score=score, mapQV=mapQV,
                nMatch=nMatch, nMismatch=nMismatch, nIns=nIns, nDel=nDel,
                qAln=qAln, alnStr=alnStr, sAln=sAln,
                strand=strand, identity=identity)
        except (AssertionError, ValueError) as e:
            errMsg = "String not recognized as a valid BLASR M5 record."
            raise ValueError(errMsg + str(e))


class BLASRReaderBase(object):

    """BLASR M4, M5 Reader Base."""

    def __init__(self, fileName, className="BLASRReaderBase"):
        self.fileName = fileName
        self.className = className
        try:
            self.infile = open(self.fileName, 'r')
        except IOError as e:
            errMsg = self.className + ": could not read file " + \
                fileName + "\n" + str(e)
            raise IOError(errMsg)

    def __iter__(self):
        raise NotImplementedError(self.className +
                                  ".__iter__ not impelemented.")

    def __enter__(self):
        return self

    def close(self):
        """Close the M5 file."""
        self.infile.close()

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()


class BLASRM4Reader(BLASRReaderBase):

    """M4 reader."""

    def __init__(self, fileName, flipQS=False):
        BLASRReaderBase.__init__(self, fileName=fileName,
                                 className="BLASRM4Reader")
        self.flipQS = flipQS

    def __iter__(self):
        try:
            for line in self.infile:
                line = line.strip()
                if len(line) == 0 or line[0] == '#' or \
                        line.lower().endswith("mapQV".lower()):
                    continue
                yield BLASRRecord.fromM4String(line=line,
                                               flipQS=self.flipQS)
        except ValueError as e:
            raise ValueError(self.className + " failed to parse "
                             + self.fileName + "\n" + str(e))


class BLASRM5Reader(BLASRReaderBase):

    """
    Reader for BLASR M5 output.
    M5 format:
    0     1       2      3    4       5     6       7      8    9
    qName qLength qStart qEnd qStrand tName tLength tStart tEnd tStrand
    10    11       12          13     14     15
    score numMatch numMismatch numIns numDel mapQV
    16          17           18
    qAlignedSeq matchPattern tAlignedSeq
    """

    def __init__(self, fileName):
        BLASRReaderBase.__init__(self, fileName=fileName,
                                 className="BLASRM5Reader")

    def __iter__(self):
        try:
            for line in self.infile:
                line = line.strip()
                if len(line) == 0 or line[0] == '#' or \
                        line.lower().endswith("tAlignedSeq".lower()):
                    continue
                yield BLASRRecord.fromM5String(line=line)
        except ValueError as e:
            raise ValueError(self.className + " failed to parse "
                             + self.fileName + "\n" + str(e))

