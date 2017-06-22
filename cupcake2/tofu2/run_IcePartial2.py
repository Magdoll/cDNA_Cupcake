__author__ = 'etseng@pacb.com'

"""
Modified ice_partial.py for ToFU2
"""

import logging
import sys
import os.path as op

from cupcake2.tofu2.ClusterOptions2 import IceOptions2, SgeOptions2

from pbcore.util.ToolRunner import PBMultiToolRunner
from pbtranscript.__init__ import get_version
from pbtranscript.ClusterOptions import SgeOptions
from pbtranscript.PBTranscriptOptions import _wrap_parser

from cupcake2.ice2.IcePartialAll2 import IceAllPartials2, add_ice_all_partials_arguments

from cupcake2.ice2.IcePartial2 import IcePartialOne2, add_ice_partial_one_arguments

from cupcake2.ice2.IcePartialSplit2 import IcePartialSplit2, add_ice_partial_split_arguments

##from pbtranscript.ice.IcePartialI import IcePartialI, \
#    add_ice_partial_i_arguments

from pbtranscript.ice.IcePartialMerge import IcePartialMerge, \
    add_ice_partial_merge_arguments


class IcePartialRunner(PBMultiToolRunner):

    """ice_partial runner, subcommands include 'all', 'one',
    'split', 'i', and 'merge'"""

    def __init__(self):
        desc = "Toolkit for assigning non-full-length reads to isoforms."
        super(IcePartialRunner, self).__init__(desc)
        subparsers = self.subParsers

        parser = subparsers.add_parser('all',
                                       description=IceAllPartials2.desc)
        add_ice_all_partials_arguments(parser)

        parser = subparsers.add_parser('one',
                                       description=IcePartialOne2.desc)
        add_ice_partial_one_arguments(parser)

        parser = subparsers.add_parser('split',
                                       description=IcePartialSplit2.desc)
        add_ice_partial_split_arguments(parser)
        #
        # parser = subparsers.add_parser('i',
        #                                description=IcePartialI.desc)
        # add_ice_partial_i_arguments(parser)
        #
        parser = subparsers.add_parser('merge',
                                        description=IcePartialMerge.desc)
        add_ice_partial_merge_arguments(parser)

    def getVersion(self):
        """Return version string."""
        return get_version()

    def run(self):
        """Execute ice_partial.py all|one|split|i|merge."""
        cmd = self.args.subCommand
        logging.info("Running {f} {cmd} v{v}.".format(f=op.basename(__file__),
                                                      cmd=cmd, v=get_version()))
        cmd_str = ""
        try:
            args = self.args
            obj = None
            if cmd in ('all', 'one'):
                # currently user NOT allowed to set full missed start/end
                # (we also set it to 30/10 bp which is more stringent than in IceIterative2,
                # which is hard set to 50/30)
                ice_opts = IceOptions2(ece_penalty=args.ece_penalty,
                           ece_min_len=args.ece_min_len,
                           max_missed_start=args.max_missed_start,
                           max_missed_end=args.max_missed_end,
                           full_missed_start=30,
                           full_missed_end=10,
                           min_match_len=50,
                           aligner_choice=args.aligner_choice)


            if cmd == "all":
                sge_opts = SgeOptions2(unique_id=args.unique_id,
                                      use_sge=args.use_sge,
                                      max_sge_jobs=args.max_sge_jobs,
                                      sge_queue=args.sge_queue,
                                      sge_env_name=args.sge_env_name,
                                      qsub_extra=args.qsub_extra)
                obj = IceAllPartials2(root_dir=args.root_dir,
                                     fasta_filenames=args.fasta_filenames.split(','),
                                     fastq_filenames=args.fastq_filenames.split(',') if args.fastq_filenames is not None else None,
                                     ref_fasta=args.ref_fasta,
                                     out_pickle=args.out_pickle,
                                     ice_opts=ice_opts,
                                     sge_opts=sge_opts,
                                     cpus=args.cpus,
                                     tmp_dir=args.tmp_dir)
            elif cmd == "one":
                # Only assign nfl reads in the given input_fasta file to isoforms
                # "one" is always run locally so no need for SGE option
                obj = IcePartialOne2(input_fasta=args.input_fasta,
                                     input_fastq=args.input_fastq,
                                     ref_fasta=args.ref_fasta,
                                     out_pickle=args.out_pickle,
                                     done_filename=args.done_filename,
                                     ice_opts=ice_opts,
                                     cpus=args.cpus,
                                     tmp_dir=args.tmp_dir)
            elif cmd == "split":
                obj = IcePartialSplit2(root_dir=args.root_dir,
                                       nfl_fa=args.nfl_fa,
                                       nfl_fq=args.nfl_fq,
                                       N=args.N)
            # elif cmd == "i":
            #     obj = IcePartialI(root_dir=args.root_dir, i=args.i,
            #                       ccs_fofn=args.ccs_fofn,
            #                       blasr_nproc=args.blasr_nproc,
            #                       tmp_dir=args.tmp_dir)
            elif cmd == "merge":
                 obj = IcePartialMerge(root_dir=args.root_dir,
                                       N=args.N)
            else:
                raise ValueError("Unknown command passed to {f}: {cmd}.".
                                 format(f=op.basename(__file__), cmd=cmd))

            cmd_str = obj.cmd_str()
            logging.info("Running CMD: {cmd_str}".format(cmd_str=cmd_str))
            obj.run()
        except:
            logging.exception("Exiting {cmd_str} with return code 1.".
                              format(cmd_str=cmd_str))
            return 1
        return 0


def main():
    """Main function."""
    runner = IcePartialRunner()
    return runner.start()

if __name__ == "__main__":
    sys.exit(main())

