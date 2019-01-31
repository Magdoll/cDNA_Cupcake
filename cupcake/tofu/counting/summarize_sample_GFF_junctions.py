__author__ = 'etseng@pacb.com'

"""
Reporter for junction summary for one or more samples.
Recommended to run before scrubbing sample GFFs prior to chaining.

Suggested process is:
1. run collapse to get GFF for each sample
2. run this report script on all sample GFFs
3. run scrubber on all sample GFFs
4. run chaining using scrubbed (cleaned) sample GFFs
"""
import os, sys
from collections import defaultdict
import cupcake.io.GFF as GFF
import numpy as np
from Bio.Seq import Seq
from Bio import SeqIO
from sklearn.cluster import Birch

def sanity_check(sample_dirs, gff_filename, genome_filename=None, junction_filename=None):
    for d in sample_dirs.itervalues():
        file = os.path.join(d, gff_filename)
        if not os.path.exists(file):
            print >> sys.stderr, "Expected GFF file {0} does not exist. Abort!".format(file)
            sys.exit(-1)

    if genome_filename is not None and not os.path.exists(genome_filename):
        print >> sys.stderr, "Genome file {0} given but does not exist. Abort!".format(genome_filename)
        sys.exit(-1)

    if junction_filename is not None and not os.path.exists(junction_filename):
        print >> sys.stderr, "Junction file {0} given but does not exist. Abort!".format(junction_filename)
        sys.exit(-1)

def read_config(filename):
    """
    SAMPLE=<name>;<path>

    must also have
    GFF_FILENAME=

    optional:
    GENOME_FILENAME=
    JUNCTION_FILENAME=
    GROUP_FILENAME=

    Everything else will be ignored (so you can re-use sample.config for chain_samples.py)
    """
    sample_dirs = {}
    sample_names = []
    gff_filename = None
    genome_filename = None
    junction_filename = None

    with open(filename) as f:
        for line in f:
            if line.startswith('tmpSAMPLE='):
                print >> sys.stderr, "Please only use SAMPLE=, not tmpSAMPLE= for junction reports!"
                sys.exit(-1)
            elif line.startswith('SAMPLE='):
                name, path = line.strip()[len('SAMPLE='):].split(';')
                if name.startswith('tmp_'):
                    print >> sys.stderr, "Sample names are not allowed to start with tmp_! Please change {0} to something else.".format(name)
                    sys.exit(-1)
                sample_dirs[name] = os.path.abspath(path)
                sample_names.append(name)
            elif line.startswith('GFF_FILENAME='):
                gff_filename = line.strip()[len('GFF_FILENAME='):]
            elif line.startswith('GENOME_FILENAME='):
                genome_filename = line.strip()[len('GENOME_FILENAME='):]
            elif line.startswith('JUNCTION_FILENAME='):
                junction_filename = line.strip()[len('JUNCTION_FILENAME='):]

    if gff_filename is None:
        raise Exception, "Expected GFF_FILENAME= but not in config file {0}! Abort.".format(filename)

    if len(sample_names) == 0:
        print >> sys.stderr, "No samples given. Exit."
        sys.exit(-1)

    return sample_dirs, sample_names, gff_filename, genome_filename, junction_filename


def read_annotation_junction_bed(junction_filename):
    """
    junction.bed is in format:

    chr, left (0-based), right (0-based), +/-

    following junction.bed format from TopHat
    http://ccb.jhu.edu/software/tophat/manual.shtml
    """
    junction = defaultdict(lambda: {}) # (chr, strand) --> (start, end)
    for line in open(junction_filename):
        chrom, left, right, strand = line.strip().split('\t')
        junction[chrom, strand][(int(left), int(right))] = 1
    return junction

def summarize_junctions(sample_dirs, sample_names, gff_filename, output_prefix, genome_d=None, junction_known=None):
    """
    1. for each sample, read all the GFF, store the junction information (both 0-based)

    """
    junc_by_chr_strand = defaultdict(lambda: defaultdict(lambda: [])) # (chr,strand) --> (donor,acceptor) --> samples it show up in (more than once possible)

    for sample_name, d in sample_dirs.iteritems():
        for r in GFF.collapseGFFReader(os.path.join(d, gff_filename)):
            n = len(r.ref_exons)
            if n == 1: continue # ignore single exon transcripts
            for i in xrange(n-1):
                donor = r.ref_exons[i].end-1 # make it 0-based
                accep = r.ref_exons[i+1].start # start is already 0-based
                junc_by_chr_strand[r.chr, r.strand][donor, accep].append(sample_name)

    # write junction report
    f1 = open(output_prefix+'.junction.bed', 'w')
    f1.write("track name=junctions description=\"{0}\" useScore=1\n".format(output_prefix))
    with open(output_prefix+'.junction_detail.txt', 'w') as f:
        f.write("chr\tleft\tright\tstrand\tnum_transcript\tnum_sample\tgenome\tannotation\tlabel\n")
        keys = junc_by_chr_strand.keys()
        keys.sort()
        for _chr, _strand in keys:
            v = junc_by_chr_strand[_chr, _strand]
            v_keys = v.keys()
            v_keys.sort()
            labels = cluster_junctions(v_keys)
            for i,(_donor, _accep) in enumerate(v_keys):
                f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t".format(_chr, _donor, _accep, _strand, len(v[_donor,_accep]), len(set(v[_donor,_accep]))))
                f1.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(_chr, _donor, _accep+1, output_prefix, len(v[_donor,_accep]), _strand))
                # if genome is given, write acceptor-donor site
                if genome_d is None or _chr not in genome_d:
                    f.write("NA\t")
                else:
                    up, down = genome_d[_chr][_donor+1:_donor+3], genome_d[_chr][_accep-2:_accep]
                    if _strand == '+':
                        f.write("{0}-{1}\t".format(str(up.seq).upper(), str(down.seq).upper()))
                    else:
                        f.write("{0}-{1}\t".format(str(down.reverse_complement().seq).upper(), str(up.reverse_complement().seq).upper()))
                # if annotation is given, check if matches with annotation
                if junction_known is None:
                    f.write("NA\n")
                else:
                    if (_chr, _strand) in junction_known and (_donor, _accep) in junction_known[_chr, _strand]:
                        f.write("Y\t")
                    else:
                        f.write("N\t")
                f.write("{c}_{s}_{lab}\n".format(c=_chr, s=_strand, lab=labels[i]))
    f1.close()

    return junc_by_chr_strand

def cluster_junctions(juncs):
    birch_model = Birch(threshold=3, n_clusters=None)
    X = np.array(juncs)
    birch_model.fit(X)

    return birch_model.labels_


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("config", help="Config filename")
    parser.add_argument("output_prefix", help="Output prefix")

    args = parser.parse_args()

    sample_dirs, sample_names, gff_filename, genome_filename, junction_filename = \
       read_config(args.config)

    sanity_check(sample_dirs, gff_filename, genome_filename, junction_filename)

    if genome_filename is not None:
        print >> sys.stderr, "Reading genome file {0}...".format(genome_filename)
        genome_d = SeqIO.to_dict(SeqIO.parse(open(genome_filename), 'fasta'))
    else:
        print >> sys.stderr, "No genome file given. Ignore."
        genome_d = None

    if junction_filename is not None:
        print >> sys.stderr, "Reading junction file {0}....".format(junction_filename)
        junction_bed = read_annotation_junction_bed(junction_filename)
    else:
        print >> sys.stderr, "No junction file given. Ignore."
        junction_bed = None


    summarize_junctions(sample_dirs, sample_names, gff_filename, args.output_prefix, genome_d, junction_bed)
