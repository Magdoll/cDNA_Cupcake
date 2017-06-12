from cupcake2.tofu2.ToFuOptions2 import add_sge_arguments, \
    add_tmp_dir_argument, add_cluster_summary_report_arguments, \
    add_ice_post_arrow_hq_lq_arguments2, \
    add_cluster_root_dir_as_positional_argument, \
    add_fofn_arguments

from cupcake2.ice2.IceArrow2 import IceArrow2
from cupcake2.ice2.IceArrowPostProcess2 import IceArrowPostProcess2
from cupcake2.ice2.__init__ import ICE_ARROW_PY

class IceArrowAll2(object):
    """
    IceArrowAll2
    """

    desc = "After assigning all non-full-length reads to unpolished " + \
           "consensus isoforms (e.g., 'run_IcePartials2.py all' is done), " + \
           "polish these isoforms by using Arrow for Sequel data " + \
           "then output high QV and low QV isoforms."

    prog = "%s all " % ICE_ARROW_PY

    def __init__(self, root_dir, subread_xml, sge_opts, ipq_opts,
                 report_fn=None, summary_fn=None, tmp_dir=None, prog_name=None):
        prog_name = prog_name if prog_name is not None else "IceArrowAll2"
        self.root_dir = root_dir
        self.subread_xml = subread_xml
        self.report_fn = report_fn
        self.summary_fn = summary_fn
        self.sge_opts = sge_opts
        self.ipq_opts = ipq_opts
        self.tmp_dir = tmp_dir

    def cmd_str(self):
        return self._cmd_str(root_dir=self.root_dir, subread_xml=self.subread_xml,
                             sge_opts=self.sge_opts,
                             ipq_opts=self.ipq_opts, report_fn=self.report_fn,
                             summary_fn=self.summary_fn, tmp_dir=self.tmp_dir)


    def _cmd_str(self, root_dir, subread_xml, sge_opts, ipq_opts,
                 report_fn, summary_fn, tmp_dir):
        """Return a cmd string. ($ICE_ARROW_PY all)."""
        cmd = self.prog + \
              "{d} ".format(d=root_dir) + \
              "--subread_xml={f} ".format(f=subread_xml)
        if tmp_dir is not None:
            cmd += "--tmp_dir={d} ".format(d=tmp_dir)
        if report_fn is not None:
            cmd += "--report={f} ".format(f=report_fn)
        if summary_fn is not None:
            cmd += "--summary={f} ".format(f=summary_fn)
        cmd += sge_opts.cmd_str(show_blasr_nproc=True, show_arrow_nproc=True)
        cmd += ipq_opts.cmd_str()
        return cmd

    def run(self):
        """Run"""
        iceq = IceArrow2(root_dir=self.root_dir, subread_xml=self.subread_xml,
                         sge_opts=self.sge_opts,
                         tmp_dir=self.tmp_dir)
        iceq.validate_inputs()
        iceq.run()

        icepq = IceArrowPostProcess2(root_dir=self.root_dir,
                                     quit_if_not_done=False,
                                     ipq_opts=self.ipq_opts)
        icepq.run()
        return 0


def add_ice_arrow_all_arguments(parser):

    arg_parser = parser #parser.arg_parser.parser
    arg_parser = add_cluster_root_dir_as_positional_argument(arg_parser)
    arg_parser = add_fofn_arguments(arg_parser, subread_xml=True)
    arg_parser = add_cluster_summary_report_arguments(arg_parser)
    arg_parser = add_ice_post_arrow_hq_lq_arguments2(arg_parser)
    arg_parser = add_sge_arguments(arg_parser, blasr_nproc=True, arrow_nproc=True)
    arg_parser = add_tmp_dir_argument(arg_parser)
    return arg_parser