__author__ = 'etseng@pacb.com'

"""
Modified ice_partial.py for ToFU2
"""

#!/usr/bin/env python
###############################################################################
# Copyright (c) 2011-2013, Pacific Biosciences of California, Inc.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
# * Neither the name of Pacific Biosciences nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
#
# NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE GRANTED BY
# THIS LICENSE.  THIS SOFTWARE IS PROVIDED BY PACIFIC BIOSCIENCES AND ITS
# CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT
# NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR
# ITS CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
# OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
# OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
# ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
###############################################################################


"""
Overview:
    pbtranscript cluster contains two main components:
    * (1) ICE (iterative clustering and error correction) to predict
      unpolished consensus isoforms.
    * (2) Polish, to use nfl reads and quiver to polish those predicted
      unpolished isoforms. Polish contains three steps:
      + (2.1) IceAllPartials (ice_partial.py all)
              Align and assign nfl reads to unpolished isoforms, and
              save results to a pickle file.
      + (2.2) IceQuiver (ice_quiver.py all)
              Call quiver to polish each isoform based on alignments
              created by mapping its associated fl and nfl reads to
              this isoform.
      + (2.3) IceQuiverPostprocess (ice_quiver.py postprocess)
              Collect and post process quiver results, and classify
              HQ/LQ isoforms.

    In order to handle subtasks by SMRTPipe instead of pbtranscript
    itself, we will refactor the polish phase including
    (2.1) (2.2) and (2.3).

    (2.1) IceAllPartials (ice_partial.py all) will be refactored to
      + (2.1.1) ice_partial.py split
                Split nfl reads into N chunks (N<=100).
      + (2.1.2) ice_partial.py i
                For each chunk i, align and assign its reads to unpolished
                isoforms and create a pickle file.
      + (2.1.3) ice_partial.py merge
                Merge pickles for all splitted chunks together to a
                big pickle.

    Thus, the task of 2.1 can either be done by running:
        ice_partial.py all ...
    or by
        first splitting the nfl reads into N chunks:
            ice_partial.py split root_dir nfl_fa N
        then assigning nfl reads in each chunk to isoforms independently:
            ice_partial.py i root_dir {i} \
                           --ccs_fofn=ccs_fofn --blasr_nproc=blasr_nproc
            , for i = 0, ..., N-1
        and finally mergeing all results to a big one:
            ice_partial.py merge root_dir N

    Hierarchy:
        pbtranscript = iceiterative

        pbtranscript --quiver = iceiterative + \
                                ice_polish.py

        ice_polish.py =  ice_make_fasta_fofn.py + \
                         ice_partial.py all + \
                         ice_quiver.py all

        ice_partial.py all = ice_partial.py split + \
                             ice_partial.py i + \
                             ice_partial.py merge

        (ice_partial.py one --> only apply ice_partial on a given input fasta)

        ice_quiver.py all = ice_quiver.py i + \
                            ice_quiver.py merge + \
                            ice_quiver.py postprocess

Alternative way to call this suite of scripts:
    python -m pbtranscript.ice_partial
"""

import logging
import sys
import os.path as op

from cupcake2.tofu2.ClusterOptions2 import IceOptions2, SgeOptions2

from pbcore.util.ToolRunner import PBMultiToolRunner
from pbtranscript.__init__ import get_version
from pbtranscript.ClusterOptions import SgeOptions
from pbtranscript.PBTranscriptOptions import _wrap_parser

from cupcake2.ice2.IceAllPartials2 import IceAllPartials2, add_ice_all_partials_arguments

from cupcake2.ice2.IcePartial2 import IcePartialOne2, add_ice_partial_one_arguments

##from pbtranscript.ice.IcePartialI import IcePartialI, \
#    add_ice_partial_i_arguments
#from pbtranscript.ice.IcePartialSplit import IcePartialSplit, \
#    add_ice_partial_split_arguments
#from pbtranscript.ice.IcePartialMerge import IcePartialMerge, \
#    add_ice_partial_merge_arguments


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

        # parser = subparsers.add_parser('split',
        #                                description=IcePartialSplit.desc)
        # add_ice_partial_split_arguments(parser)
        #
        # parser = subparsers.add_parser('i',
        #                                description=IcePartialI.desc)
        # add_ice_partial_i_arguments(parser)
        #
        # parser = subparsers.add_parser('merge',
        #                                description=IcePartialMerge.desc)
        # add_ice_partial_merge_arguments(parser)

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

            ice_opts = IceOptions2(ece_penalty=args.ece_penalty,
                       ece_min_len=args.ece_min_len,
                       max_missed_start=args.max_missed_start,
                       max_missed_end=args.max_missed_end,
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
                                     fastq_filenames=args.fastq_filenames.split(','),
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
            # elif cmd == "split":
            #     obj = IcePartialSplit(root_dir=args.root_dir,
            #                           nfl_fa=args.nfl_fa,
            #                           N=args.N)
            # elif cmd == "i":
            #     obj = IcePartialI(root_dir=args.root_dir, i=args.i,
            #                       ccs_fofn=args.ccs_fofn,
            #                       blasr_nproc=args.blasr_nproc,
            #                       tmp_dir=args.tmp_dir)
            # elif cmd == "merge":
            #     obj = IcePartialMerge(root_dir=args.root_dir,
            #                           N=args.N)
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

