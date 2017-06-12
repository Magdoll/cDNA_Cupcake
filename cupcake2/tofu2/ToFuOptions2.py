import random
import argparse

from pbcommand.cli.core import get_default_argparser_with_base_opts
from pbcommand.models import FileTypes, SymbolTypes, ResourceTypes, get_pbparser, PbParser
from pbtranscript.Utils import validate_fofn



class BaseConstants(object):
    TOOL_ID = "pbtranscript.tasks.main"
    DRIVER_EXE = "pbtranscript --resolved-tool-contract"
    PARSER_DESC = ""

    # these are all referenced in 2.3 pipeline scripts
    MIN_SEQ_LEN_ID = "pbtranscript.task_options.min_seq_len"
    MIN_SEQ_LEN_DEFAULT = 50
    #MIN_SCORE_ID = "pbtranscript.task_options.min_score"
    MIN_SCORE_DEFAULT = 10
    REQUIRE_POLYA_ID = "pbtranscript.task_options.require_polya"
    REQUIRE_POLYA_DEFAULT = True
    PRIMER_SEQUENCES_ID = "pbtranscript.task_options.primer_sequences"
    PRIMER_SEQUENCES_DEFAULT = ""

    HQ_ARROW_CIDS_PER_FILE = 100  # we put 100 cids in one .bam file to run Arrow

    HQ_ARROW_MIN_ACCURACY_ID = "pbtranscript.task_options.hq_arrow_min_accuracy"
    HQ_ARROW_MIN_ACCURACY_DEFAULT = 0.99
    HQ_ARROW_MIN_ACCURACY_DESC = "Minimum allowed Arrow accuracy to classify an isoform " + \
                                  "as hiqh-quality (default: %s)" % HQ_ARROW_MIN_ACCURACY_DEFAULT

    QV_TRIM_FIVEPRIME_ID = "pbtranscript.task_options.qv_trim_5p"
    QV_TRIM_FIVEPRIME_DEFAULT = 100
    QV_TRIM_FIVEPRIME_DESC = "Ignore QV of n bases in the 5' end " + \
            "(default %s)." % QV_TRIM_FIVEPRIME_DEFAULT
    QV_TRIM_THREEPRIME_ID = "pbtranscript.task_options.qv_trim_3p"
    QV_TRIM_THREEPRIME_DEFAULT = 30
    QV_TRIM_THREEPRIME_DESC = "Ignore QV of n bases in the 3' end " + \
            "(default %s)." % QV_TRIM_THREEPRIME_DEFAULT
    SAMPLE_NAME_ID = "pbtranscript.task_options.sample_name"
    SAMPLE_NAME_DEFAULT = ""

    USE_FINER_QV_ID = "pbtranscript.task_options.use_finer_qv"
    USE_FINER_QV_DEFAULT = False


def add_sge_arguments(arg_parser, blasr_nproc=False, arrow_nproc=False, gcon_nproc=False):
    """Add Sge arguments as a group to parser, return parser."""
    sge_group = arg_parser.add_argument_group("SGE environment arguments")

    sge_group.add_argument("--use_sge",
                           dest="use_sge",
                           default=False,
                           action="store_true",
                           help="Use SGE computing cluster (default: False)")

    sge_group.add_argument("--max_sge_jobs",
                           type=int,
                           dest="max_sge_jobs",
                           default=40,
                           action="store",
                           help="The maximum number of jobs that will " +
                                "be submitted to SGE concurrently. (default: 40)")

    sge_group.add_argument("--unique_id",
                           type=int,
                           dest="unique_id",
                           action="store",
                           default=random.randint(1, 100000000),
                           help="Unique ID for submitting SGE jobs. (default: randomly generated)")
    if blasr_nproc is True:
        sge_group.add_argument("--blasr_nproc",
                               type=int,
                               dest="blasr_nproc",
                               action="store",
                               default=12,
                               help="Number of cores for each BLASR. (default: 12, ignored if using DALIGNER)")
    if arrow_nproc is True:
        sge_group.add_argument("--arrow_nproc",
                               dest="arrow_nproc",
                               type=int,
                               default=8,
                               help="Number of CPUs each Arrow job uses. (default: 8)")
    if gcon_nproc is True:
        sge_group.add_argument("--gcon_nproc",
                               dest="gcon_nproc",
                               type=int,
                               default=4,
                               help="Number of CPUs for each PBDagcon job. (default: 4)")

    sge_group.add_argument("--sge_env_name",
                           type=str,
                           dest="sge_env_name",
                           default="smp",
                           action="store",
                           help="SGE parallel environment, e.g, smp, orte (default: smp).")

    sge_group.add_argument("--sge_queue",
                           type=str,
                           dest="sge_queue",
                           default=None,
                           action="store",
                           help="SGE queue to submit jobs. (default: None)")

    sge_group.add_argument("--qsub_extra",
                           type=str,
                           dest="qsub_extra",
                           default='',
                           action="store",
                           help="Additional qsub parameters")

    return arg_parser


def add_ice_arguments(arg_parser, tc_parser=None):
    """Add Ice options as a group to parser, return parser"""
    ice_group = arg_parser.add_argument_group("ICE arguments")

    # ice_group.add_argument("--quiver",
    #                        dest="quiver",
    #                        default=False,
    #                        action="store_true",
    #                        help="Call quiver to polish consensus isoforms " +
    #                             "using non-full-length non-chimeric CCS "+
    #                             "reads.")

    ice_group.add_argument("--aligner_choice", default='daligner', choices=['daligner', 'blasr'], help="Aligner choice (default: daligner)")
    ice_group.add_argument("--targeted_isoseq",
                           dest="targeted_isoseq",
                           default=False,
                           action="store_true",
                           help="Input data is from targeted Iso-Seq. Automatically make parameter adjustments.")

    ice_group.add_argument("--ece_penalty", dest="ece_penalty", type=int, default=1)
    ice_group.add_argument("--ece_min_len", dest="ece_min_len", type=int, default=20)

    ice_group.add_argument("--max_missed_start", type=int, default=400, help="Max missed start (default: 400 bp)")
    ice_group.add_argument("--max_missed_end", type=int, default=50, help="Max missed end (default: 50 bp)")

    # number of flnc reads per split
    ice_group.add_argument("--flnc_reads_per_split",
                           type=int,
                           action="store",
                           default=20000,
                           help=argparse.SUPPRESS)
    # number of nfl reads per split
    # ice_group.add_argument("--nfl_reads_per_split",
    #                        type=int,
    #                        action="store",
    #                        default=30000,
    #                        help=argparse.SUPPRESS)
    #
    # desc = "Use finer classes of QV information from CCS input instead of "+\
    #        "a single QV from FASTQ.  This option is slower and consumes "+\
    #        "more memory."
    # ice_group.add_argument("--use_finer_qv",
    #                        dest="use_finer_qv",
    #                        default=False,
    #                        action="store_true",
    #                        help=desc)


    return arg_parser



def add_fofn_arguments(arg_parser, subread_xml=False, tool_contract_parser=None):
    helpstr = "A FOFN of bam or dataset xml (e.g., input.fofn|bam|subreadset.xml), " + \
              "which contain quality values of raw reads and subreads"
    if subread_xml is True:
        arg_parser.add_argument("--subread_xml",
                                dest="subread_xml",
                                type=validate_fofn,
                                default=None,
                                action="store",
                                help=helpstr)
        if tool_contract_parser is not None:
            tool_contract_parser.add_input_file_type(FileTypes.DS_SUBREADS,
                "subreads_fofn", "SubreadSet", helpstr)

#    helpstr = "A FOFN of trimmed subreads fasta (e.g. input.fasta.fofn)"
#    if fasta_fofn is True:
#        arg_parser.add_argument("--fasta_fofn",
#                                dest="fasta_fofn",
#                                type=validate_fofn,
#                                default=None,
#                                help=helpstr)
    return arg_parser


def add_tmp_dir_argument(parser):
    """Add an argument for specifying tempoary directory, which
    will be created for storing intermediate files of isoseq
    clusters.
    """
    helpstr = "Directory to store temporary files." + \
              "(default, write to root_dir/tmp.)."
    parser.add_argument("--tmp_dir", default=None, type=str,
                        dest="tmp_dir", help=helpstr)
    return parser


def add_partial_argument(parser):
    helpstr = "Aligner choices and match parameters"
    parser.add_argument("--aligner_choice", type=str, choices=['blasr', 'daligner'], default='daligner', help="Aligner choices (default: daligner)")
    parser.add_argument("--cpus", type=int, default=4, help="Number of CPUs aligner uses (default: 4)")
    parser.add_argument("--max_missed_start", type=int, default=400, help="Max missed start (NOTE: currently ignored)")
    parser.add_argument("--max_missed_end", type=int, default=50, help="Max missed end (NOTE: currently ignored)")
    parser.add_argument("--ece_penalty", type=int, default=1, help="ECE penalty (default: 1)")
    parser.add_argument("--ece_min_len", type=int, default=20, help="ECE min length (default: 20bp)")
    return parser


def add_nfl_fa_argument(arg_parser, positional=False, required=False,
        tool_contract_parser=None):
    """
    Add nfl_fa or --nfl_fa, can be positional  or non-positional,
    required or not required.  (For the tool contract interface, it is
    always required, since we will be guaranteed to have this file internally.)
    """
    # note however that the nfl transcripts may not be used by the cluster
    # command if polishing mode is not run - they will however be handled by
    # the separately run ice_partial task
    helpstr = "Input non-full-length reads in FASTA or ContigSet format, used for " + \
              "polishing consensus isoforms, e.g., isoseq_nfl.fasta"
    if positional is True:
        arg_parser.add_argument("nfl_fa", type=str, help=helpstr)
    else:
        assert(required is True or required is False)
        arg_parser.add_argument("--nfl_fa", type=str, dest="nfl_fa",
                                required=required, help=helpstr)

    if tool_contract_parser is not None:
        tool_contract_parser.add_input_file_type(FileTypes.DS_CONTIG,
            "nfl_fa", "FASTA or ContigSet file", helpstr)
    return arg_parser


def add_nfl_fq_argument(arg_parser, positional=False, required=False,
        tool_contract_parser=None):
    helpstr = "Input non-full-length reads in FASTQ or ContigSet format, used for " + \
              "polishing consensus isoforms, e.g., isoseq_nfl.fastq"
    if positional is True:
        arg_parser.add_argument("nfl_fq", type=str, help=helpstr)
    else:
        assert(required is True or required is False)
        arg_parser.add_argument("--nfl_fq", type=str, dest="nfl_fq",
                                required=required, help=helpstr)

    if tool_contract_parser is not None:
        tool_contract_parser.add_input_file_type(FileTypes.DS_CONTIG,
            "nfl_fq", "FASTQ or ContigSet file", helpstr)
    return arg_parser


def add_cluster_root_dir_as_positional_argument(arg_parser):
    """Add root_dir as root output directory for `cluster`."""
    helpstr = "An directory to store temporary and output cluster files" + \
              " (e.g., in SA3.2, tasks/pbtranscript.tasks.separate_flnc-0/0to1kb_part0/cluster_out; in SA2.3, data/clusterOutDir/)"
    arg_parser.add_argument("root_dir", help=helpstr, type=str)
    return arg_parser

def add_cluster_summary_report_arguments(parser):
    """Add a report and a summary to summarize isoseq cluster."""
    helpstr = "CSV file to output cluster info " + \
              "(default: *.cluster_report.csv)"
    #p1 = parser.tool_contract_parser
    #p2 = parser.arg_parser.parser
    p2 = parser

    helpstr = "TXT file to output cluster summary " + \
              "(default: *.cluster_summary.txt)"
    p2.add_argument("--summary", default=None, type=str,
                        dest="summary_fn", help=helpstr)
    #p1.add_output_file_type(FileTypes.JSON, "json_summary",
    #    name="Transcript Clustering Report",
    #    description="JSON summary",
    #    default_name="summary")

    # FIXME make this a REPORT instead?
    #p1.add_output_file_type(FileTypes.CSV, "cluster_report",
    #    name="Clustering Results",
    #    description="Clustering results for each CCS read",
    #    default_name="cluster_report")
    p2.add_argument("--report", default=None, type=str,
                    dest="report_fn", help=helpstr)

    p2.add_argument("--pickle_fn", default=None, type=str,
        dest="pickle_fn")
    return parser




def add_ice_post_arrow_hq_lq_arguments2(parser):
    """Add quiver QV threshold to mark an isoform as high-quality or low-quality."""
    # if isinstance(parser, PbParser):
    #     #parser = _wrap_parser(parser)
    #     arg_parser = parser.arg_parser.parser
    #     tcp = parser.tool_contract_parser
    #     tcp.add_float(BaseConstants.HQ_ARROW_MIN_ACCURACY_ID, "hq_arrow_min_accuracy",
    #                   default=BaseConstants.HQ_ARROW_MIN_ACCURACY_DEFAULT,
    #                   name="Minimum Quiver|Arrow Accuracy", description=BaseConstants.HQ_ARROW_MIN_ACCURACY_DESC)
    #     tcp.add_int(BaseConstants.QV_TRIM_FIVEPRIME_ID, "qv_trim_5",
    #                 default=BaseConstants.QV_TRIM_FIVEPRIME_DEFAULT,
    #                 name="Trim QVs 5'", description=BaseConstants.QV_TRIM_FIVEPRIME_DESC)
    #     tcp.add_int(BaseConstants.QV_TRIM_THREEPRIME_ID, "qv_trim_3",
    #                 default=BaseConstants.QV_TRIM_THREEPRIME_DEFAULT,
    #                 name="Trim QVs 3'", description=BaseConstants.QV_TRIM_THREEPRIME_DESC)
    # else:
    #     assert isinstance(parser, argparse.ArgumentParser)
    #     arg_parser = parser

    arg_parser = parser
    icq_gp = arg_parser.add_argument_group("IceArrow High QV/Low QV arguments")
    icq_gp.add_argument("--hq_arrow_min_accuracy",
                        type=float,
                        default=BaseConstants.HQ_ARROW_MIN_ACCURACY_DEFAULT,
                        dest="hq_arrow_min_accuracy",
                        help=BaseConstants.HQ_ARROW_MIN_ACCURACY_DESC)
    icq_gp.add_argument("--qv_trim_5",
                        type=int,
                        default=BaseConstants.QV_TRIM_FIVEPRIME_DEFAULT,
                        dest="qv_trim_5",
                        help=BaseConstants.QV_TRIM_FIVEPRIME_DESC)
    icq_gp.add_argument("--qv_trim_3",
                        type=int,
                        default=BaseConstants.QV_TRIM_THREEPRIME_DEFAULT,
                        dest="qv_trim_3",
                        help=BaseConstants.QV_TRIM_THREEPRIME_DESC)


    icq_gp = arg_parser.add_argument_group("IceArrow2 HQ/LQ IO arguments")
    icq_gp.add_argument("--hq_isoforms_fa",
                        default=None,
                        type=str,
                        dest="hq_isoforms_fa",
                        help="Arrow polished, high quality isoforms " +
                        "in FASTA, default: root_dir/output/all_arrowed_hq.fasta")

    icq_gp.add_argument("--hq_isoforms_fq",
                        default=None,
                        type=str,
                        dest="hq_isoforms_fq",
                        help="Arrow polished, high quality isoforms " +
                        "in FASTQ, default: root_dir/output/all_arrowed_hq.fastq")

    icq_gp.add_argument("--lq_isoforms_fa",
                        default=None,
                        type=str,
                        dest="lq_isoforms_fa",
                        help="Arrow polished, low quality isoforms " +
                        "in FASTA, default: root_dir/output/all_arrowed_lq.fasta")

    icq_gp.add_argument("--lq_isoforms_fq",
                        default=None,
                        type=str,
                        dest="lq_isoforms_fq",
                        help="Arrow polished, low quality isoforms " +
                        "in FASTQ, default: root_dir/output/all_arrowed_lq.fastq")
    return parser
