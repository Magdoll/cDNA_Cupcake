import os, sys
from csv import DictReader
from pbtranscript.Utils import real_upath

def generate_batch_cmds(csv_filename, dirname, cmd_filename, nfl_filename, tucked_filename, subread_xml, cpus):
    cmd_f = open(cmd_filename, 'w')
    for r in DictReader(open(csv_filename), delimiter=','):
        cid = r['cluster']
        d2 = os.path.join(dirname, cid)
        if not os.path.exists(d2):
            print >> sys.stderr, "Directory {0} does not exist! Abort!".format(d2)
            sys.exit(-1)

        cmd_f.write("cd {0}\n".format(real_upath(d2)))
        cmd_f.write("fa2fq.py isoseq_flnc.fasta\n")
        cmd_f.write("run_IceInit2.py isoseq_flnc.fasta init.uc.pickle --aligner_choice=blasr --cpus={c}\n".format(c=cpus))
        cmd_f.write("run_IceIterative2.py isoseq_flnc.fasta isoseq_flnc.fastq isoseq_flnc.fasta . " \
                    "--init_uc_pickle=init.uc.pickle --aligner_choice=blasr " + \
                    "--blasr_nproc {c} --gcon_nproc {c2}\n".format(c=cpus, c2=min(cpus, 4)))
        cmd_f.write("run_IcePartial2.py all {nfl},{tucked} ".format(nfl=nfl_filename, tucked=tucked_filename) + \
                    "output/final.consensus.fasta nfl.pickle " + \
                    "--root_dir . --aligner_choice=blasr --cpus={c}\n".format(c=cpus))
        cmd_f.write("run_IceArrow2.py all --subread_xml {s} ".format(s=subread_xml) + \
                    "--blasr_nproc {c} --arrow_nproc {c} .\n".format(c=cpus))

    cmd_f.close()






if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser("Generate batch commands for running Ice2->Partial2->Arrow2 for each preCluster output bin")
    parser.add_argument("precluster_csv", help="Cluster CSV file (ex: preCluster.cluster_info.csv)")
    parser.add_argument("precluster_dir", help="preCluster out directory (ex: preCluster_out/)")
    parser.add_argument("nfl_filename", help="nFL filename (ex: isoseq_nfl.fasta)")
    parser.add_argument("tucked_filename", help="tucked filename (ex: preCluster_out.tucked.fasta)")
    parser.add_argument("subread_xml", help="Subread XML")
    parser.add_argument("--cpus", default=20, type=int, help="Number of CPUs (default: 20)")
    parser.add_argument("--cmd_filename", default='cmds', help="Output command filename (default: cmds)")
    args = parser.parse_args()

    generate_batch_cmds(args.precluster_csv, real_upath(args.precluster_dir),
                        args.cmd_filename,
                        real_upath(args.nfl_filename),
                        real_upath(args.tucked_filename),
                        real_upath(args.subread_xml),
                        args.cpus)