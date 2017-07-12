
import os, sys
from csv import DictReader
from Bio import SeqIO
from pbtranscript.Utils import real_upath
from cupcake2.io.SeqSplitter import splitFaFq

"""
USAGE:

generate_batch_cmd_for_preCluster_out.py [csv] [dir] [cmd]

where [csv] is the preCluster_out.cluster_info.csv and
      [dir] is the preCluster_out/ directory with all the bins.

[csv] should contain <cluster> and <size>.

Generates a series of commands that can be run either locally or split for qsub.
Only runs the IceInit2 --> IceIterative2 step! Not the IcePartial2->IceArrow2 step!
"""


def fa2fq(input):
    """
    Copy of fa2fq from sequence/
    """
    try:
        assert input.lower().endswith('.fasta') or input.lower().endswith('.fa')
    except AssertionError:
        print >> sys.stderr, "Input {0} does not end with .fasta or .fa! Abort".format(input)
        sys.exit(-1)
    output = input[:input.rfind('.')] + '.fastq'

    f = open(output, 'w')
    for r in SeqIO.parse(open(input),'fasta'):
        r.letter_annotations['phred_quality'] = [60]*len(r.seq)
        SeqIO.write(r, f, 'fastq')
    f.close()

    print >> sys.stderr, "Output written to", f.name
    return f.name


def preprocess_flnc_split_if_necessary(dirname, flnc_size, flnc_split=20000):
    """
    Go into a precluster bin and look at isoseq_flnc.fasta
    If it is bigger than <flnc_split> sequences, than split and return the splits
    """
    if flnc_size <= flnc_split:
        fa_files = [os.path.join(dirname, 'isoseq_flnc.fasta')]
    else:
        fa_files = splitFaFq(os.path.join(dirname, 'isoseq_flnc.fasta'),
                  reads_per_split=flnc_split,
                  out_dir=dirname,
                  out_format="input_split{0:03d}.fasta",
                  is_fq=False,
                  reads_in_first_split=flnc_split/2)
    # convert all fasta to fastq
    fq_files = [fa2fq(x) for x in fa_files]
    return [os.path.basename(x) for x in fa_files], [os.path.basename(x) for x in fq_files]


def generate_batch_cmds(csv_filename, dirname, cmd_filename, cpus):
    #, nfl_filename, tucked_filename, subread_xml, cpus):
    cmd_f = open(cmd_filename, 'w')
    for r in DictReader(open(csv_filename), delimiter=','):
        cid = r['cluster']
        d2 = os.path.join(dirname, cid)
        if not os.path.exists(d2):
            print >> sys.stderr, "Directory {0} does not exist! Abort!".format(d2)
            sys.exit(-1)

        cmd_f.write("cd {0}\n".format(real_upath(d2)))

        fa_files, fq_files = preprocess_flnc_split_if_necessary(d2, int(r['size']), flnc_split=20000)

        cmd_f.write("run_IceInit2.py {fa} init.uc.pickle --aligner_choice=blasr --cpus={c}\n".format(c=cpus, fa=fa_files[0]))
        cmd_f.write("run_IceIterative2.py {fas} {fqs} isoseq_flnc.fasta . ".format(fas=",".join(fa_files), fqs=",".join(fq_files)) + \
                    "--init_uc_pickle=init.uc.pickle --aligner_choice=blasr " + \
                    "--blasr_nproc {c} --gcon_nproc {c2}\n".format(c=cpus, c2=min(cpus, 4)))
#        cmd_f.write("run_IcePartial2.py all {nfl},{tucked} ".format(nfl=nfl_filename, tucked=tucked_filename) + \
#                    "output/final.consensus.fasta nfl.pickle " + \
#                    "--root_dir . --aligner_choice=blasr --cpus={c}\n".format(c=cpus))
#        cmd_f.write("run_IceArrow2.py all --subread_xml {s} ".format(s=subread_xml) + \
#                    "--blasr_nproc {c} --arrow_nproc {c} .\n".format(c=cpus))

    cmd_f.close()






if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser("Generate batch commands for running IceInit2->IceIterative2 for each preCluster output bin")
    parser.add_argument("precluster_csv", help="Cluster CSV file (ex: preCluster.cluster_info.csv)")
    parser.add_argument("precluster_dir", help="preCluster out directory (ex: preCluster_out/)")
    #parser.add_argument("nfl_filename", help="nFL filename (ex: isoseq_nfl.fasta)")
    #parser.add_argument("tucked_filename", help="tucked filename (ex: preCluster_out.tucked.fasta)")
    #parser.add_argument("subread_xml", help="Subread XML")
    parser.add_argument("--cpus", default=20, type=int, help="Number of CPUs (default: 20)")
    parser.add_argument("--cmd_filename", default='cmds', help="Output command filename (default: cmds)")
    args = parser.parse_args()

    generate_batch_cmds(args.precluster_csv, real_upath(args.precluster_dir),
                        args.cmd_filename,
                        #real_upath(args.nfl_filename),
                        #real_upath(args.tucked_filename),
                        #real_upath(args.subread_xml),
                        args.cpus)