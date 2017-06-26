__author__ = 'etseng@pacb.com'

#!/usr/bin/env python

import logging
import sys
import os.path as op
from pbcore.util.ToolRunner import PBMultiToolRunner
from pbtranscript.__init__ import get_version


from cupcake2.tofu2.ClusterOptions2 import SgeOptions2, IceArrowHQLQOptions2

from cupcake2.ice2.IceArrowAll2 import IceArrowAll2, add_ice_arrow_all_arguments

from cupcake2.ice2.IceArrowPostProcess2 import IceArrowPostProcess2, add_ice_arrow_postprocess_arguments

# from pbtranscript.ice.IceQuiverI import IceQuiverI, \
#     add_ice_quiver_i_arguments
# from pbtranscript.ice.IceQuiverMerge import IceQuiverMerge, \
#     add_ice_quiver_merge_arguments
# from pbtranscript.ice.IceQuiverPostprocess import \
#     IceQuiverPostprocess, add_ice_quiver_postprocess_arguments


class IceArrowRunner2(PBMultiToolRunner):

    """ice_quiver runner, subcommands include 'all', 'i', 'merge' and
    'postprocess' """

    def __init__(self):
        desc = "Toolkit for assigning non-full-length reads to isoforms."
        super(IceArrowRunner2, self).__init__(desc)
        subparsers = self.subParsers

        parser = subparsers.add_parser('all',
                                       description=IceArrowAll2.desc)
        add_ice_arrow_all_arguments(parser)

        # parser = subparsers.add_parser('i',
        #                                description=IceQuiverI.desc)
        # add_ice_quiver_i_arguments(parser)
        #
        # parser = subparsers.add_parser('merge',
        #                                description=IceQuiverMerge.desc)
        # add_ice_quiver_merge_arguments(parser)
        #
        parser = subparsers.add_parser('postprocess',
                                        description=IceArrowPostProcess2.desc)
        add_ice_arrow_postprocess_arguments(parser)

    def getVersion(self):
        """Return version string."""
        return get_version()

    def run(self):
        """Execute ice_quiver.py all|i|merge|postprocess."""
        cmd = self.args.subCommand
        logging.info("Running {f} {cmd} v{v}.".format(f=op.basename(__file__),
                                                      cmd=cmd, v=get_version()))
        cmd_str = ""
        try:
            args = self.args
            obj = None
            if cmd == "all":
                sge_opts = SgeOptions2(unique_id=args.unique_id,
                                       use_sge=args.use_sge,
                                       max_sge_jobs=args.max_sge_jobs,
                                       blasr_nproc=args.blasr_nproc,
                                       arrow_nproc=args.arrow_nproc,
                                       sge_env_name=args.sge_env_name,
                                       sge_queue=args.sge_queue,
                                       qsub_extra=args.qsub_extra)
                ipq_opts = IceArrowHQLQOptions2(
                    hq_isoforms_fa=args.hq_isoforms_fa,
                    hq_isoforms_fq=args.hq_isoforms_fq,
                    lq_isoforms_fa=args.lq_isoforms_fa,
                    lq_isoforms_fq=args.lq_isoforms_fq,
                    qv_trim_5=args.qv_trim_5,
                    qv_trim_3=args.qv_trim_3,
                    hq_arrow_min_accuracy=args.hq_arrow_min_accuracy)
                obj = IceArrowAll2(root_dir=args.root_dir,
                                   subread_xml=args.subread_xml,
                                   sge_opts=sge_opts,
                                   ipq_opts=ipq_opts,
                                   tmp_dir=args.tmp_dir)
            # elif cmd == "i":
            #     sge_opts = SgeOptions(unique_id=args.unique_id,
            #                           use_sge=args.use_sge,
            #                           max_sge_jobs=args.max_sge_jobs,
            #                           blasr_nproc=args.blasr_nproc,
            #                           arrow_nproc=args.arrow_nproc)
            #     obj = IceQuiverI(root_dir=args.root_dir, i=args.i, N=args.N,
            #                      bas_fofn=args.bas_fofn,
            #                      fasta_fofn=None,
            #                      sge_opts=sge_opts,
            #                      tmp_dir=args.tmp_dir)
            #elif cmd == "merge":
            #    obj = IceQuiverMerge(root_dir=args.root_dir, N=args.N)
            elif cmd == "postprocess":
                ipq_opts = IceArrowHQLQOptions2(
                    hq_isoforms_fa=args.hq_isoforms_fa,
                    hq_isoforms_fq=args.hq_isoforms_fq,
                    lq_isoforms_fa=args.lq_isoforms_fa,
                    lq_isoforms_fq=args.lq_isoforms_fq,
                    qv_trim_5=args.qv_trim_5,
                    qv_trim_3=args.qv_trim_3,
                    hq_arrow_min_accuracy=args.hq_arrow_min_accuracy,
                    hq_min_full_length_reads=args.hq_min_full_length_reads)
                obj = IceArrowPostProcess2(root_dir=args.root_dir,
                                           ipq_opts=ipq_opts,
                                           quit_if_not_done=args.quit_if_not_done,
                                           summary_fn=args.summary_fn,
                                           report_fn=args.report_fn)
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
    runner = IceArrowRunner2()
    return runner.start()

if __name__ == "__main__":
    sys.exit(main())
