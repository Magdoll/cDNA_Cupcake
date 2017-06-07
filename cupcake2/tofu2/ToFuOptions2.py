import random

from pbcommand.cli.core import get_default_argparser_with_base_opts
from pbcommand.models import FileTypes, SymbolTypes, ResourceTypes, get_pbparser, PbParser
from pbtranscript.Utils import validate_fofn

def add_sge_arguments(arg_parser, blasr_nproc=False, quiver_nproc=False, gcon_nproc=False):
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
    if quiver_nproc is True:
        sge_group.add_argument("--quiver_nproc",
                               dest="quiver_nproc",
                               type=int,
                               default=8,
                               help="Number of CPUs each quiver job uses. (default: 8)")
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



def add_fofn_arguments(arg_parser, bas_fofn=False, fasta_fofn=False,
                       tool_contract_parser=None):
    """Add bas_fofn, fasta_fofn arguments."""
    helpstr = "A FOFN of bam or dataset xml (e.g., input.fofn|bam|subreadset.xml), " + \
              "which contain quality values of raw reads and subreads"
    if bas_fofn is True:
        arg_parser.add_argument("--subread_xml",
                                dest="subreads_xml",
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
    parser.add_argument("--max_missed_end", type=int, default=100, help="Max missed end (NOTE: currently ignored)")
    parser.add_argument("--ece_penalty", type=int, default=1, help="ECE penalty (default: 1)")
    parser.add_argument("--ece_min_len", type=int, default=20, help="ECE min length (default: 20bp)")
    return parser