__author__ = 'etseng@pacb.com'

from Bio.Seq import Seq
from Bio.SeqIO import SeqRecord

class LazyFastqReader:
    """
    Like LazyFastaReader except works with fastq!
    """
    def __init__(self, fastq_filename, seqid_extraction=lambda x: x):
        self.f = open(fastq_filename)
        self.d = {}

        while 1:
            line = self.f.readline()
            if len(line) == 0: break
            assert line.startswith('@')
            id = line.strip()[1:].split(None, 1)[0] # the header MUST be just 1 line
            id2 = seqid_extraction(id)
            if id2 in self.d:
                raise Exception("Duplicate id {0} (or {1})!!".format(id, id2))
            self.d[id2] = self.f.tell()
            self.f.readline() # seq
            self.f.readline() # +
            self.f.readline() # quality

    def __getitem__(self, k):
        if k not in self.d:
            raise Exception("key {0} not in dictionary!".format(k))
        self.f.seek(self.d[k])

        sequence = self.f.readline().strip()
        assert self.f.readline().startswith('+')
        qualstr = self.f.readline().strip()
        return SeqRecord(seq=Seq(sequence), id=k, \
                         letter_annotations={'phred_quality':[ord(x)-33 for x in qualstr]})

    def __setitem__(self, key, value):
        self.d[key] = value

    def keys(self):
        return self.d.keys()


class LazyFastaReader:
    """
    NOTE: this version works BioPython

    This is meant to substitute for the Bio.SeqIO.to_dict method since some fasta files
    are too big to fit entirely to memory. The only requirement is that every id line
    begins with the symbol >. It is ok for the sequences to stretch multiple lines.
    The sequences, when read, are returned as FastaRecord objects.

    Example:
        r = FastaReader('output/test.fna')
        r['6C_49273_NC_008578/2259031-2259297'] ==> this shows the FastaRecord
    """
    def __init__(self, fasta_filename, seqid_extraction=lambda x: x):
        self.f = open(fasta_filename)
        self.d = {}

        while 1:
            line = self.f.readline()
            if len(line) == 0: break
            if line.startswith('>'):
                id = line.strip()[1:].split(None, 1)[0] # the header MUST be just 1 line
                id2 = seqid_extraction(id)
                if id2 in self.d:
                    raise Exception("Duplicate id {0} (or {1})!!".format(id, id2))
                self.d[id2] = self.f.tell()

    def __getitem__(self, k):
        if k not in self.d:
            raise Exception("key {0} not in dictionary!".format(k))
        self.f.seek(self.d[k])
        content = ''
        for line in self.f:
            if line.startswith('>'):
                break
            content += line.strip()
        return SeqRecord(seq=Seq(content), id=k)

    def __setitem__(self, key, value):
        self.d[key] = value

    def keys(self):
        return self.d.keys()
