
import logging
import sys
import os.path as op
from pbcore.util.ToolRunner import PBToolRunner
from pbtranscript.__init__ import get_version

from cupcake2.tofu2.ToFuOptions2 import add_cluster_root_dir_as_positional_argument
from cupcake2.ice2.__init__ import ICE_ARROW_PY

from pbtranscript.ice.IceQuiver import IceQuiver
from pbtranscript.Utils import cat_files, nfs_exists



def add_ice_quiver_merge_arguments(parser):
    """Add IceQuiverMerge arguments."""
    parser = add_cluster_root_dir_as_positional_argument(parser)

    helpstr = "Number of workload chunks."
    parser.add_argument("N", help=helpstr, type=int)
    return parser


class IceQuiverMerge(object):

    """Merge all quiver polished isoforms done by N IceQuiverI jobs."""

    desc = "Unpolished isoforms were divided into N chunks and " + \
           "polished using Quiver or Arrow separately. Now collect all " + \
           "submitted polishing jobs from " + \
           "root_dir/log/submitted_quiver_jobs.{i}of{N}.txt " + \
           "(i=0,...,N-1) to root_dir/log/submitted_quiver_jobs.txt"

    prog = "%s merge" % ICE_ARROW_PY

    def __init__(self, root_dir, N):
        self.root_dir = root_dir
        self.N = int(N)

    def getVersion(self):
        """Return version string."""
        return get_version()

    def cmd_str(self):
        """Return a cmd string of IceQuiverMerge ($ICE_ARROW_PY merge)"""
        return self._cmd_str(root_dir=self.root_dir, N=self.N)

    def _cmd_str(self, root_dir, N):
        """Return a cmd string of IceQuiverMerge ($ICE_ARROW_PY merge)."""
        cmd = " ".join([str(x) for x in (self.prog, root_dir, N)])
        return cmd

    def run(self):
        """Run"""
        iceq = IceQuiver(root_dir=self.root_dir, bas_fofn=None,
                         fasta_fofn=None, sge_opts=None,
                         prog_name="ice_quiver_merge")

        iceq.add_log(self.cmd_str())
        iceq.add_log("root_dir: {d}.".format(d=self.root_dir))
        iceq.add_log("Total number of chunks: N = {N}.".format(N=self.N))

        src = [iceq.submitted_quiver_jobs_log_of_chunk_i(i=i, num_chunks=self.N)
               for i in range(0, self.N)]
        for f in src:
            if not nfs_exists(f):
                raise IOError("Log {f} ".format(f=f) +
                              "of submitted quiver jobs does not exist.")

        dst = iceq.submitted_quiver_jobs_log

        iceq.add_log("Collecting submitted quiver jobs from:\n{src}\nto {dst}.".
                     format(src="\n".join(src), dst=dst))

        cat_files(src=src, dst=dst)

        iceq.close_log()


class IceQuiverMergeRunner(PBToolRunner):

    """IceQuiverMerge Runner"""

    def __init__(self):
        PBToolRunner.__init__(self, IceQuiverMerge.desc)
        add_ice_quiver_merge_arguments(self.parser)

    def getVersion(self):
        """Return version string."""
        return get_version()

    def run(self):
        """Run"""
        logging.info("Running {f} v{v}.".format(f=op.basename(__file__),
                                                v=get_version()))
        cmd_str = ""
        try:
            args = self.args
            iceqm = IceQuiverMerge(root_dir=args.root_dir, N=args.N)
            cmd_str = iceqm.cmd_str()
            iceqm.run()
        except:
            logging.exception("Exiting {cmd_str} with return code 1.".
                              format(cmd_str=cmd_str))
            return 1
        return 0


def main():
    """Main function."""
    runner = IceQuiverMergeRunner()
    return runner.start()

if __name__ == "__main__":
    sys.exit(main())