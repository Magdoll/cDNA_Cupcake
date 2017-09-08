__author__ = 'etseng@pacb.com'

"""
File Handling for Heng Li's minimap
https://github.com/lh3/minimap

Output format:
0. Query sequence name
1. Query sequence length
2. Query start coordinate (0-based)
3. Query end coordinate (0-based) # Liz: will convert to 1-based
4. + if query and target on the same strand; - if opposite Target sequence name
5. Target sequence name
6. Target sequence length
7. Target start coordinate on the original strand
8. Target end coordinate on the original strand
9. Number of matching bases in the mapping
10. Number bases, including gaps, in the mapping
11. Mapping quality (0-255 with 255 missing unavailable)
"""

from cupcake.io.BLASRRecord import BLASRReaderBase

class MiniRecord:
    def __init__(self, qID, qLength, qStart, qEnd, strand,
                 sID, sLength, sStart, sEnd,
                 nMatch=None, nBase=None):
        self.qID = qID
        self.qLength = qLength
        self.qStart = qStart # 0-based
        self.qEnd = qEnd     # 1-based
        self.strand = strand
        self.sID = sID
        self.sLength = sLength
        self.sStart = sStart
        self.sEnd = sEnd
        self.nMatch = nMatch
        self.nBase = nBase

        #self.status = None # unassigned

    def __str__(self):
        return """
        qID: {0}
        sID: {1}
        strand: {2}
        qStart-qEnd: {3}-{4}
        sStart-sEnd: {5}-{6}
        """.format(self.qID, self.sID, self.strand, self.qStart, self.qEnd, self.sStart, self.sEnd)


    @classmethod
    def fromPAF(cls, line):
        """
        parsing Heng Li's Pairwise mapping format (PAF), Table 2 in paper
        """
        raw = line.strip().split('\t')
        try:
            assert raw[4] in ('+', '-')
            return MiniRecord(qID=raw[0], qLength=int(raw[1]), qStart=int(raw[2]), qEnd=int(raw[3])+1,
                              strand=raw[4],
                              sID=raw[5], sLength=int(raw[6]), sStart=int(raw[7]), sEnd=int(raw[8])+1,
                              nMatch=int(raw[9]), nBase=int(raw[10]))
        except:
            raise ValueError, "String not recognized as a valid PAF record: {0}".format(line)


    def characterize(self, max_missed_5_len, max_missed_5_ratio, max_missed_3_len, max_missed_3_ratio, \
                     max_contained_len_diff, max_contained_len_diff_ratio):
        """
        q_contained: query is fully contained in subject
        s_contained: subject is fully contained in query
        <--- all above is conditioned on not having too much unmapped --->
        otherwise it is 'partial'
        """
        q_5_mapped = self.qStart <= min(max_missed_5_len, max_missed_5_ratio*self.qLength)
        q_3_mapped = self.qLength-self.qEnd <= min(max_missed_3_len, max_missed_3_ratio*self.qLength)
        s_5_mapped = self.sStart <= min(max_missed_5_len, max_missed_5_ratio*self.sLength)
        s_3_mapped = self.sLength-self.sEnd <= min(max_missed_3_len, max_missed_3_ratio*self.sLength)
        q_fully_mapped = q_5_mapped and q_3_mapped
        s_fully_mapped = s_5_mapped and s_3_mapped

        q_map_len = self.qEnd-self.qStart
        s_map_len = self.sEnd-self.sStart
        q_s_len_match = abs(q_map_len-s_map_len) <= \
                        min(max_contained_len_diff, max_contained_len_diff_ratio*max(q_map_len,s_map_len))


        if not q_fully_mapped and not s_fully_mapped: # neither is fully mapped, so partial!
            return "partial"

        if q_fully_mapped:
            if s_fully_mapped:
                if q_s_len_match:
                    return "match"
                else: # probably diff isoforms
                    return "partial"
            else:
                # q fully mapped, s is not, check if query is fully contained in subject
                # this means subject is not allowed to have too much unmapped 3' but 5' is ok
                if s_3_mapped:
                    if q_s_len_match:
                        return "q_contained"
                    else:
                        return "partial"
                else: return "partial"
        else:
            if s_fully_mapped:
                # s is fully mapped, q is not, check if subject is fully contained in query
                # this means query is not allowed to have too much unmapped 3' but 5' is ok
                if q_3_mapped:
                    if q_s_len_match:
                        return "s_contained"
                    else:
                        return "partial"
                else: return "partial"
            else: # neither is fully mapped, so partial!
                return "partial"



class MiniReader(BLASRReaderBase):
    def __init__(self, fileName):
        BLASRReaderBase.__init__(self, fileName=fileName,
                                 className="MiniReader")

    def __iter__(self):
        return self

    def next(self):
        cur = self.infile.tell()
        line = self.infile.readline().strip()
        if self.infile.tell() == cur:
            raise StopIteration, "EOF reached!!"

        return MiniRecord.fromPAF(line)

