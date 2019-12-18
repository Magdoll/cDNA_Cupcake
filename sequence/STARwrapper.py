#!/usr/bin/env python
__author__ = 'etseng@pacb.com'
"""
Wrapper for running STARlong.
Parameters are pre-set according to:


"""
import os, sys, subprocess
import tempfile
import shutil

CMD_STARlong = '/home/UNIXHOME/etseng/software_downloads/STAR-2.5.3a/bin/Linux_x86_64/STAR --runMode alignReads --outSAMattributes NH HI NM MD --readNameSeparator space --outFilterMultimapScoreRange 1 --outFilterMismatchNmax 2000 --scoreGapNoncan -1  --scoreGapGCAG -4 --scoreGapATAC -8 --scoreDelOpen -1 --scoreDelBase -1 --scoreInsOpen -1 --scoreInsBase -1 --alignEndsType Local --seedSearchStartLmax 50 --seedPerReadNmax 100000 --seedPerWindowNmax 1000 --alignTranscriptsPerReadNmax 100000 --alignTranscriptsPerWindowNmax 10000'
CMD_STAR2_format = CMD_STARlong + " --twopassMode None --runThreadN {c} --genomeDir {d} --readFilesIn {i}"

def run_STAR(args):
    tmp_dir= tempfile.mkdtemp(prefix='STARtmp')

    in_fasta = os.path.abspath(args.in_fasta)
    out_sam = os.path.abspath(args.out_sam)
    cur_dir = os.getcwd()
    os.chdir(tmp_dir)
    cmd = CMD_STAR2_format.format(c=args.cpus, d=args.genome_dir, i=in_fasta)
    if subprocess.check_call(cmd, shell=True)!=0:
        print("ERROR RUNNING CMD:", cmd, file=sys.stderr)
        sys.exit(-1)

    os.chdir(cur_dir)
    shutil.move(os.path.join(tmp_dir, 'Aligned.out.sam'), out_sam)
    shutil.rmtree(tmp_dir)

if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser("Wrapper for running STARlong.")
    parser.add_argument("genome_dir")
    parser.add_argument("in_fasta")
    parser.add_argument("out_sam")
    parser.add_argument("--cpus", type=int, default=10, help="Number of threads (default: 10)")


    args = parser.parse_args()
    run_STAR(args)
