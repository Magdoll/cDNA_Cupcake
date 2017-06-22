__author__ = 'etseng@pacb.com'

"""
ice_opts = IceOptions2(aligner_choice='daligner', max_missed_start=400, max_missed_end=50)
icec = pp.IceIterative(fasta_filename='isoseq_flnc.fasta', fasta_filenames_to_add=[], all_fasta_filename='isoseq_flnc.fasta', ccs_fofn=None, root_dir='.', ice_opts=ice_opts,
sge_opts=sge_opts, uc=a, probQV=None, fastq_filename='isoseq_flnc.fastq', output_pickle_file='test.pickle', tmp_dir='./tmp')
icec.run()
"""

import random
from cPickle import load

from cupcake2.tofu2.ToFuOptions2 import add_ice_arguments, add_sge_arguments,\
    add_cluster_root_dir_as_positional_argument
from cupcake2.tofu2.ClusterOptions2 import SgeOptions2, IceOptions2

from cupcake2.ice2.IceIterative2 import IceIterative2
from cupcake2.tofu2.run_IceInit2 import run_IceInit2

def run_IceIterative2(fasta_filenames, fastq_filenames, all_fasta_filename, ice_opts, sge_opts, root_dir, init_uc_pickle=None):
    """

    """
    if init_uc_pickle is not None:
        uc = load(open(init_uc_pickle))
    else:
        uc = run_IceInit2(readsFa=fasta_filenames[0],
                        out_pickle=init_uc_pickle,
                        ice_opts=ice_opts,
                        sge_opts=sge_opts,
                        readsFq=fastq_filenames[0])

    icec = IceIterative2(fasta_filename=fasta_filenames[0],
                         fastq_filename=fastq_filenames[0],
                         fasta_filenames_to_add=fasta_filenames[1:],
                         fastq_filenames_to_add=fastq_filenames[1:],
                         all_fasta_filename=all_fasta_filename,
                         root_dir=root_dir,
                         ice_opts=ice_opts,
                         sge_opts=sge_opts,
                         uc=uc,
                         probQV=None,
                         refs=None,
                         d=None,
                         is_FL=True,
                         qv_prob_threshold=0.03,
                         output_pickle_file=None,
                         tmp_dir=None)
    icec.run()


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser("run IceIterative2.")

    parser.add_argument("fasta_filenames")
    parser.add_argument("fastq_filenames")
    parser.add_argument("all_fasta_filename")

    parser = add_cluster_root_dir_as_positional_argument(parser)
    parser.add_argument("--init_uc_pickle", default=None, type=str, help="Init2 pickle (ex: init.uc.pickle)")

    parser = add_ice_arguments(parser)
    parser = add_sge_arguments(parser, blasr_nproc=True, gcon_nproc=True, arrow_nproc=False)


    args = parser.parse_args()

    # currently, the user is NOT allowed to set the definition for full missed start/end
    # (i.e. how much of the query must be mapped to consider it a hit)
    ice_opts = IceOptions2(ece_penalty=args.ece_penalty,
                           ece_min_len=args.ece_min_len,
                           max_missed_start=args.max_missed_start,
                           max_missed_end=args.max_missed_end,
                           full_missed_start=50,
                           full_missed_end=30,
                           aligner_choice=args.aligner_choice,
                           )

    unique_id = random.randint(1, 1000000)
    sge_opts = SgeOptions2(unique_id,
                          use_sge=args.use_sge,
                          max_sge_jobs=args.max_sge_jobs,
                          blasr_nproc=args.blasr_nproc,
                          gcon_nproc=args.gcon_nproc,
                          sge_queue=args.sge_queue,
                          sge_env_name=args.sge_env_name,
                          qsub_extra=args.qsub_extra)

    run_IceIterative2(args.fasta_filenames.split(','),
                      args.fastq_filenames.split(','),
                      args.all_fasta_filename,
                      ice_opts,
                      sge_opts,
                      args.root_dir,
                      args.init_uc_pickle)
