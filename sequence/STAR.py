#!/usr/bin/env python

"""
Readers for STAR (spliced aligner) format reading
"""

class STARJunctionRecord:
    """
    column 1: chromosome
    column 2: first base of the intron (1-based) --> store as 0-based
    column 3: last base of the intron (1-based)
    column 4: strand (0: undefined, 1: +, 2: -)
    column 5: intron motif: 0: non-canonical; 1: GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5: AT/AC, 6: GT/AT
    column 6: 0: unannotated, 1: annotated (only if splice junctions database is used) 10
    column 7: number of uniquely mapping reads crossing the junction
    column 8: number of multi-mapping reads crossing the junction
    column 9: maximum spliced alignment overhang
    """
    strand_dict = {0: 'NA', 1: '+', 2: '-'}
    motif_dict = ['non_canonical', 'GTAG', 'CTAC', 'GCAG', 'CTGC', 'ATAC', 'GTAT']

    def __init__(self, chrom, start, end, strand, motif, is_annotated, unique_count, multi_count, overhang):
        self.chrom = chrom
        self.start = start  # store as 0-based start
        self.end = end      # store as 1-based end
        self.strand = strand
        self.motif = motif
        self.is_annotated = is_annotated
        self.unique_count = unique_count
        self.multi_count = multi_count
        self.overhang = overhang

    def __str__(self):
        return """
        Junction: {c}:{s}-{e} ({t})
        Motif: {m}
        Annotated: {a}
        Counts: {u} (unique), {n} (multi)
        """.format(c=self.chrom, s=self.start+1, e=self.end, t=self.strand,
                   m=self.motif, a=self.is_annotated, u=self.unique_count,
                   n=self.multi_count)

    @staticmethod
    def process_line(line):
        raw = line.strip().split()
        if len(raw)!=9:
            raise Exception("Expected 9 columns for STAR junction file! Got {0} instead!".format(len(raw)))

        chrom = raw[0]
        start = int(raw[1])
        end = int(raw[2])
        strand = STARJunctionRecord.strand_dict[int(raw[3])]


        return STARJunctionRecord(chrom=raw[0],
                                  start=int(raw[1])-1,
                                  end=int(raw[2]),
                                  strand=STARJunctionRecord.strand_dict[int(raw[3])],
                                  motif=STARJunctionRecord.motif_dict[int(raw[4])],
                                  is_annotated=True if int(raw[5])==1 else False,
                                  unique_count=int(raw[6]),
                                  multi_count=int(raw[7]),
                                  overhang=int(raw[8]))


class STARJunctionReader:
    def __init__(self, filename,):
        self.filename = filename
        self.f = open(filename)

    def __iter__(self):
        return self

    def __next__(self):
        line = self.f.readline().strip()
        if len(line) == 0:
            raise StopIteration
        return STARJunctionRecord.process_line(line)
