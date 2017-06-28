__author__ = 'etseng@pacb.com'

from cPickle import dump

from pbtranscript.ice.IceUtils import set_probqv_from_model, set_probqv_from_fq
from cupcake2.ice2.IceInit2 import IceInit2
from cupcake2.tofu2.ClusterOptions2 import IceOptions2, SgeOptions2

def run_IceInit2(readsFa, out_pickle, ice_opts, sge_opts, readsFq=None):
    if readsFq is None:
        _probqv, msg = set_probqv_from_model()
    else:
        _probqv, msg = set_probqv_from_fq(readsFq)

    qver_get_func=  _probqv.get_smoothed
    qvmean_get_func= _probqv.get_mean

    i = IceInit2(readsFa, qver_get_func, qvmean_get_func, ice_opts, sge_opts)

    with open(out_pickle, 'w') as f:
        dump(i.uc, f)

    return i.uc


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser("Run IceInit2: initial clique finding on full-length reads.")
    parser.add_argument("flnc_fa", help="FLNC fasta filename")
    parser.add_argument("output_pickle", help="Output pickle filename")
    parser.add_argument("--flnc_fq", default=None, help="FLNC fastq fileame (to provide QV information)")
    parser.add_argument("--ece_penalty", default=1, type=int, help="ECE penalty (default: 1)")
    parser.add_argument("--ece_min_len", default=20, type=int, help="ECE min len (default: 20)")
    parser.add_argument("--max_missed_start", default=100, type=int, help="Max missed 5'(default: 100 bp)" )
    parser.add_argument("--max_missed_end", default=30, type=int, help="Max missed 3' (default: 30 bp)")
    parser.add_argument("--aligner_choice", default='daligner', choices=['daligner', 'blasr'], help="Aligner (default: daligner)")
    parser.add_argument("--cpus", default=4, help="Number of CPUs aligner uses (default: 4)")

    args = parser.parse_args()

    ice_opts = IceOptions2(ece_penalty=args.ece_penalty,
                           ece_min_len=args.ece_min_len,
                           max_missed_start=args.max_missed_start,
                           max_missed_end=args.max_missed_end,
                           full_missed_start=30,
                           full_missed_end=20,
                           aligner_choice=args.aligner_choice,
                           )

    sge_opts = SgeOptions2(unique_id=123, blasr_nproc=args.cpus)

    # set min_match_len to (low_cDNA_size-max_missed_start-max_missed_end), rounded to lowest 100 bp
    ice_opts.detect_cDNA_size(args.flnc_fa)
    x = ice_opts.low_cDNA_size - ice_opts.max_missed_start - ice_opts.max_missed_end
    ice_opts.min_match_len = max(50, (x/100)*100)

    run_IceInit2(args.flnc_fa, args.output_pickle, ice_opts, sge_opts, args.flnc_fq)