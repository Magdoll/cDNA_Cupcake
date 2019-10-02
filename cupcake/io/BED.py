__author__ = 'etseng@pacb.com'

"""
Misc APIs for reading BED format

UCSC BED format:
https://genome.ucsc.edu/FAQ/FAQformat#format1

SimpleBED --- chr, 0-based start, 1-based end
"""

class SimpleBED(object):
    def __init__(self, chrom, start, end, name=None):
        self.chr = chrom
        self.start = start
        self.end = end
        self.name = name

    def __str__(self):
        return "{c}:{s}-{e} (name:{n})".format(c=self.chr, s=self.start, e=self.end, n=self.name)


class SimpleBEDReader:
    def __init__(self, filename, start_base=0, end_base=1):
        self.filename = filename
        self.f = open(filename)
        
        if start_base!=0 and start_base!=1:
            raise Exception("start_base can only be 0 or 1!")
        if end_base!=0 and end_base!=1:
            raise Exception("end_base can only be 0 or 1!")

        self.start_base = start_base
        self.end_base = end_base

    def __iter__(self):
        return self

    def __next__(self):
        return self.read()

    def read(self):
        cur = self.f.tell()
        line = self.f.readline()
        if self.f.tell() == cur:
            raise StopIteration

        raw = line.strip().split('\t')
        if len(raw) >= 4: name=raw[3]
        else: name=None
        return SimpleBED(raw[0], int(raw[1])-self.start_base, int(raw[2])+(1-self.end_base), name)
