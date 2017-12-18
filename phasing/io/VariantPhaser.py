__author__ = 'etseng@pacb.com'

import pdb
from collections import defaultdict, namedtuple
from csv import DictReader
import vcf
from Bio.Seq import Seq
from Bio import SeqIO
from cupcake.io.BioReaders import GMAPSAMReader
from coordinate_mapper import get_base_to_base_mapping_from_sam


__VCF_EXAMPLE__ = \
"""
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT
20      1       .       G       A,T     .       PASS    AF=0.5;DB       GT
"""

def type_fa_or_fq(file):
    file = file.upper()
    if file.endswith('.FA') or file.endswith('.FASTA'): return 'fasta'
    else: return 'fastq'


class VariantPhaser(object):
    def __init__(self, vc):
        """
        :param vc: MPileUPVariant instance.
        """
        self.vc = vc
        self.min_var_pos = min(vc.variant)  # mininum 0-based position of a called variant
        self.max_var_pos = max(vc.variant)  # maximum 0-based position of a called variant
        self.accepted_vars_by_pos = {} # 0-based pos --> list of accepted, (NOT strand sense) base
        self.count_of_vars_by_pos = {} # 0-based pos --> (NOT strand sense, but ref-based) base --> count
        self.accepted_pos = [] # sorted list of variant positions (0-based, ref)

        # process vc.variant which is
        # dict of 0-based pos --> desc list of (base, count)
        # ex: {1565: [('a', 49), ('g', 36)]}
        # lower case means at pos 1565, we expect - strand mapping and
        # seq base is 'T' on the sense strand
        # this converts to self.accepted_vars_by_pos[1565] = ['A', 'G']
        # later, when we are matchin back to transcript seq, need to watch for strand!
        for pos, vars in vc.variant.iteritems():
            self.accepted_vars_by_pos[pos] = [_base.upper() for _base,_count in vars]
            self.count_of_vars_by_pos[pos] = dict((_base.upper(), _count) for _base,_count in vars)

        self.accepted_pos = self.accepted_vars_by_pos.keys()
        self.accepted_pos.sort()

        self.haplotypes = Haplotypes(self.accepted_pos, self.vc.ref_base, self.count_of_vars_by_pos)
        self.seq_hap_info = {} # haplotype assignment, key: (CCS) seqid, value: haplotype index


    def phase_variant(self, sam_filename, input_fa_or_fq, output_prefix, partial_ok=False):
        """
        :param sam_filename: CCS SAM filename. Can be unsorted.
        :param input_fa_or_fq: Input CCS fasta/fastq filename.
        :param output_prefix: Output prefix. Writes to xxx.log.
        :param partial_ok: default False. if True, (CCS) reads don't need to cover all SNP positions.

        For each alignment:
        1. discard if did not map to the strand expected
        2. discard if did not map to the full range of variants (unless <partial_ok> is True)
        3. discard if at var positions have non-called bases (outliers)
        """
        f_log = open(output_prefix+'.log', 'w')

        seq_dict = SeqIO.to_dict(SeqIO.parse(open(input_fa_or_fq), type_fa_or_fq(input_fa_or_fq)))
        for r in GMAPSAMReader(sam_filename, True, query_len_dict=dict((k, len(seq_dict[k].seq)) for k in seq_dict)):
            if r.sID == '*':
                f_log.write("Ignore {0} because: unmapped.\n".format(r.qID))
                continue
            if r.flag.strand != self.vc.expected_strand:
                f_log.write("Ignore {0} because: strand is {1}.\n".format(r.qID, r.flag.strand))
                continue # ignore
            if not partial_ok and (r.sStart > self.min_var_pos or r.sEnd < self.max_var_pos):
                f_log.write("Ignore {0} because: aln too short, from {1}-{2}.\n".format(r.qID, r.sStart+1, r.sEnd))
                continue

            i, msg = self.match_haplotype(r, str(seq_dict[r.qID].seq).upper(), partial_ok)
            if i is None: # read is rejected for reason listed in <msg>
                f_log.write("Ignore {0} because: {1}.\n".format(r.qID, msg))
                continue
            else:
                f_log.write("{0} phased: haplotype {1}={2}\n".format(r.qID, i, self.haplotypes[i]))
                print "{0} has haplotype {1}:{2}".format(r.qID, i, self.haplotypes[i])
                self.seq_hap_info[r.qID] = i


    def match_haplotype(self, r, s, partial_ok=False):
        """
        Match an alignment record to existing haplotypes or create a new one.
        Helper function for self.phase_variant()
        :param r: CCS alignment (SAM record)
        :param s: CCS sequence (in strand), must be plain str and every base is upper case
        :param partial_ok: default False. if True, (CCS) reads don't need to cover all SNP positions.

        :return: (haplotype_index, msg) or (None, msg) if variants don't match w/ called SNPs
        """
        assert type(s) is str and str.isupper(s)
        assert r.flag.strand == self.vc.expected_strand
        # m: mapping of 0-based seq --> 0-based ref position
        # rev_map: mapping of 0-based ref position --> 0-based seq
        m = get_base_to_base_mapping_from_sam(r.segments, r.cigar, r.qStart, r.qEnd, r.flag.strand)
        ref_m = dict((v,k) for k,v in m.iteritems())

        # go through each variant
        # <hap> to represent the concatenated string of all variant positions for this seq
        # ex: if there are three var positions, a hap would be "ATG" or "A?G" (if partial_ok is True), etc.
        hap = ''
        impute_later = False
        for ref_pos in self.accepted_pos:
            if ref_pos not in ref_m:
                if partial_ok: # read does not cover one of the SNP positions, so use "?"
                    hap += "?"
                else:
                    return None, "Does not have base at ref_pos {0}.\n".format(ref_pos)
            else:
                base = s[ref_m[ref_pos]]
                if self.vc.expected_strand == '-': # must convert the base to the rev comp
                    base = str(Seq(base).reverse_complement()).upper()
                if base in self.accepted_vars_by_pos[ref_pos]:
                    hap += base
                else: # contains a base at a variant position that is not called. Try to impute.
                    hap += base
                    impute_later = True

        if all(b=='?' for b in hap):
            return None, "Does not cover any variant base."

        if impute_later:
            impute_i = self.haplotypes.impute_haplotype(hap, min_score=3)
            if impute_i is None:
                return None, "Seq {0} contained non-called variant. Impute failed.\n".format(hap)
            else:
                return impute_i, "IMPUTED"
        return self.haplotypes.match_or_add_haplotype(hap_string=hap)



def phase_isoforms(read_stat_filename, seqids, phaser):
    """
    :param read_stat_filename: the .read_stat file that has columns <id> and <pbid>, where <id> is CCS id and <pbid> is PB.X.Y
    :param seqids: CCS IDs that were used to create the haplotypes.
    :param phaser: VariantPhaser object that contains the haplotype and seqid->haplotype information.

    :return: list of (isoform, dict of haplotype count), ex: {'PB.45.1': {0:10, 1:20}}
             which means PB.45.1 has haplotype 0 supported by 10 CCS reads and hap 1 supported by 20 CCS reads.

    *NOTE* currently uses FL CCS reads only (even if the SNPs may have been called by FL+nFL CCS SAM)
    """
    result = {} # dict of (isoform, dict of haplotype_index --> CCS count supporting it
    # from read stat, gather which isoforms have which (CCS) seq members.
    isoforms = defaultdict(lambda: []) # key: PB.X.Y, value: list of seqid members
    for r in DictReader(open(read_stat_filename), delimiter='\t'):
        if r['id'] in seqids and r['is_fl']=='Y':
            isoforms[r['pbid']].append(r['id'])

    # for each isoform, look at the CCS membership to know which haplotypes are expressed
    for _iso, _seqids in isoforms.iteritems():
        tally = defaultdict(lambda: 0) # haplotype index --> count (of CCS)
        for seqid in _seqids:
            if seqid in phaser.seq_hap_info: # some CCS (seqids) may not have been used by the phaser, so account for that
                tally[phaser.seq_hap_info[seqid]] += 1
        if len(tally) > 0:
            result[_iso] = dict(tally)
    return result


class Haplotypes(object):
    """
    Storing haplotypes for a loci.

    self.haplotype[i] is the i-th haplotype.
    if N = len(self.haplotype[i]), then there are N variants along the loci.
    self.hap_var_positions[j] means that the j-th variant corressponds to (0-based) position on the ref genome.
    """
    def __init__(self, var_positions, ref_at_pos, count_of_vars_by_pos):
        """
        :param var_positions: sorted list of (0-based) variant positions
        :param ref_at_pos: dict of (0-based) variant position --> ref base at this position
        :param count_of_vars_by_pos: 0-based pos --> (NOT strand sense, but ref-based) base --> count
        """
        self.haplotypes = [] # haplotypes, where haplotypes[i] is the i-th distinct haplotype of all var concat
        self.hap_var_positions = var_positions
        self.ref_at_pos = ref_at_pos # dict of (0-based) pos --> ref base
        self.alt_at_pos = None # init: None, later: dict of (0-based) pos --> unique list of alt bases
        self.count_of_vars_by_pos = count_of_vars_by_pos
        self.haplotype_vcf_index = None # init: None, later: dict of (hap index) --> (0-based) var pos --> phase (0 for ref, 1+ for alt)

        # sanity check: all variant positions must be present
        self.sanity_check()

    def __getitem__(self, ith):
        """
        Returns the <i>-th haplotype
        """
        return self.haplotypes[ith]

    def __str__(self):
        return """
        var positions: {pp}
        haplotypes: \n{h}
        """.format(pp=",".join(map(str,self.hap_var_positions)),
                   h="\n".join(self.haplotypes))

    def sanity_check(self):
        """
        Sanity check the following:
        -- variant positions are properly recorded and concordant
        -- alt bases are truly alt and unique
        -- all haplotypes are the same length
        """
        for pos in self.hap_var_positions:
            assert pos in self.ref_at_pos

        if self.alt_at_pos is not None:
            for pos in self.alt_at_pos:
                # ref base must not be in alt
                assert self.ref_at_pos[pos] not in self.alt_at_pos[pos]
                # alt bases must be unique
                assert len(self.alt_at_pos[pos]) == len(set(self.alt_at_pos[pos]))

        if len(self.haplotypes) >= 1:
            n = len(self.haplotypes[0])
            assert n == len(self.hap_var_positions)
            for hap_str in self.haplotypes[1:]:
                assert len(hap_str) == n


    def match_or_add_haplotype(self, hap_string):
        """
        If <hap_string> is an existing haplotype, return the index.
        Otherwise, add to known haplotypes and return the new index.

        :return: <index>, "FOUND" or "NEW"
        """
        if hap_string in self.haplotypes:
            i = self.haplotypes.index(hap_string)
            return i, "FOUND"
        else:
            i = len(self.haplotypes)
            self.haplotypes.append(hap_string)
            return i, "NEW"

    def impute_haplotype(self, hap_string, min_score):
        """
        :param hap_string: a hap string with '?'s.
        :param min_sim: minimum similarity with existing haplotype to accept assignment
        :return: <index> of an existing haplotype, or None if not sufficiently matched

        Impute haplotype and only return a match if:
        (a) score (similarity) is >= min_score
        (b) the matching score for the best one is higher than the second best match
        """
        sim_tuple = namedtuple('sim_tuple', 'index score')
        sims = [] # list of sim_tuple
        hap_str_len = len(hap_string)
        for i in xrange(len(self.haplotypes)):
            # Liz note: currently NOT checking whether existing haplotypes have '?'. I'm assuming no '?'.
            score = sum((hap_string[k]==self.haplotypes[i][k]) for k in xrange(hap_str_len))
            if score > 0:
                sims.append(sim_tuple(index=i, score=score))
        if len(sims) == 0:
            return None
        sims.sort(key=lambda x: x.score, reverse=True)
        if sims[0].score >= min_score and (len(sims)==1 or sims[0].score > sims[1].score):
            return sims[0].index
        else:
            return None

    def get_haplotype_vcf_assignment(self):
        """
        Must be called before self.write_haplotype_to_vcf()
        This is preparing for writing out VCF. We need to know, for each variant position,
        the ref base (already filled in self.ref_at_pos) and the alt bases (self.alt_at_pos).
        For each haplotype in (self.haplotype), we need to know the whether the i-th variant is the
        ref (index 0), or some alt base (index 1 and onwards).

        Propagates two variables:

        self.haplotype_vcf_index: hap index --> pos --> phase index (0 for ref, 1+ for alt)
        self.alt_at_pos: dict of <0-based pos> --> alt bases (not is not ref) at this position
        """
        self.haplotype_vcf_index = [{} for i in xrange(len(self.haplotypes))]
        self.alt_at_pos = {}

        # what happens in the case of partial phasing
        # ex: self.haplotypes[0] = "A?G", this means when it comes to the second pos, pos2,
        # in the VCF we would want to write out .|. for diploid, . for haploid, etc
        # so let's set self.haplotype_vcf_index[0][pos2] = '.' to indicate that

        for i,pos in enumerate(self.hap_var_positions):
            ref = self.ref_at_pos[pos]
            # need to go through the haplotype bases, if ref is already represented, then don't put it in alt
            self.alt_at_pos[pos] = []
            for hap_i, hap_str in enumerate(self.haplotypes):
                base = hap_str[i]
                if base=='?': # means this haplotype does not cover this position!
                    self.haplotype_vcf_index[hap_i][pos] = '.'
                elif base==ref: # is the ref base
                    self.haplotype_vcf_index[hap_i][pos] = 0
                else: # is an alt base, see if it's already there
                    if base in self.alt_at_pos[pos]:
                        j = self.alt_at_pos[pos].index(base)
                        self.haplotype_vcf_index[hap_i][pos] = j + 1 # always +1, buz alt starts at 1 (0 is ref)
                    else:
                        j = len(self.alt_at_pos[pos])
                        self.alt_at_pos[pos].append(base)
                        self.haplotype_vcf_index[hap_i][pos] = j + 1 # always +1, buz alt starts at 1 (0 is ref)
            # in the case where partial_ok=False, it's possible some alt are never presented by a haplotype
            # we must check that all variants are presented here
            for _base in self.count_of_vars_by_pos[pos]:
                if (_base not in self.ref_at_pos[pos]) and (_base not in self.alt_at_pos[pos]):
                    self.alt_at_pos[pos].append(_base)


    def write_haplotype_to_vcf(self, fake_genome_mapping_filename, isoform_tally, output_prefix):
        """
        The following functions must first be called first:
        -- self.get_haplotype_vcf_assignment
        """
        if self.haplotype_vcf_index is None or self.alt_at_pos is None:
            raise Exception, "Must call self.get_haplotype_vcf_assignment() first!"

        self.sanity_check()

        name_isoforms = isoform_tally.keys()
        name_isoforms.sort()

        # write a fake VCF example so we can read the headers in
        with open('template.vcf', 'w') as f:
            f.write(__VCF_EXAMPLE__)
        reader = vcf.VCFReader(open('template.vcf'))
        reader.samples = name_isoforms
        f_vcf = vcf.Writer(open(output_prefix+'.vcf', 'w'), reader)


        # human readable text:
        # first line: assoc VCF filename
        # second line: haplotype, list of sorted isoforms
        # third line onwards: haplotype and assoc count
        f_human = open(output_prefix+'.human_readable.txt', 'w')
        f_human.write("Associated VCF file: {0}.vcf\n".format(output_prefix))
        f_human.write("haplotype\t{samples}\n".format(samples="\t".join(name_isoforms)))
        for hap_index,hap_str in enumerate(self.haplotypes):
            f_human.write(hap_str)
            for _iso in name_isoforms:
                if hap_index in isoform_tally[_iso]:
                    f_human.write("\t{0}".format(isoform_tally[_iso][hap_index]))
                else:
                    f_human.write("\t0")
            f_human.write('\n')
        f_human.close()


        # read fake genome mapping file
        fake_map = {} # 0-based position on fake --> (chr, 0-based ref position)
        with open(fake_genome_mapping_filename) as f:
            for line in f:
                fake_pos, ref_chr, ref_pos = line.strip().split(',')
                fake_map[int(fake_pos)] = (ref_chr, int(ref_pos))


        # for each position, write out the ref and alt bases
        # then fill in for each isoform (aka "sample"):
        #  if this isoform only shows one allele, then it's just that allele (0 for ref, 1+ otherwise)
        #  if this isoform shows 2+ allele, then the first allele is indicated by self.haplotypes[0]
        for i,pos in enumerate(self.hap_var_positions):
            ref_chr, ref_pos = fake_map[pos]
            total_count = sum(self.count_of_vars_by_pos[pos].itervalues())
            alt_freq = ["{0:.2f}".format(self.count_of_vars_by_pos[pos][b]*1./total_count) for b in self.alt_at_pos[pos]]
            rec = vcf.model._Record(CHROM=ref_chr,
                              POS=ref_pos+1,
                              ID='.',
                              REF=self.ref_at_pos[pos],
                              ALT=[vcf.model._Substitution(b) for b in self.alt_at_pos[pos]],
                              QUAL='.',
                              FILTER='PASS',
                              INFO={'AF':alt_freq, 'DP':total_count},
                              FORMAT="GT:HQ",
                              sample_indexes=None)
            samp_ft = vcf.model.make_calldata_tuple(['GT', 'HQ'])
            rec.samples = []
            for _iso in name_isoforms:
                # isoform_tally[_iso] is a dict of haplotype index --> count
                # the index for thos base at this pos would thus be haplotype_vcf_index[hap_index][i]
                # we always need to show the phases in haplotype index order sorted
                hap_indices = isoform_tally[_iso].keys()
                hap_indices.sort()
                genotype = "|".join(str(self.haplotype_vcf_index[hap_index][pos]) for hap_index in hap_indices)
                counts = ",".join(str(isoform_tally[_iso][hap_index]) for hap_index in hap_indices)
                rec.samples.append(vcf.model._Call(rec, _iso, samp_ft(*[genotype, counts])))
            f_vcf.write_record(rec)
        f_vcf.close()