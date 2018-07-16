__author__ = 'etseng@pacb.com'

#!/usr/bin/env python
import os, sys
import subprocess

def run_cmd(cmd):
    if subprocess.check_call(cmd, shell=True)!=0:
        raise Exception, "Error cmd:", cmd


input = sys.argv[1] # ex: test.sam

if not input.endswith('.sam'):
    print >> sys.stderr, "Only accepts files ending in .sam. Abort!"
    sys.exit(-1)

prefix = input[:-4]

run_cmd("samtools view -bS {0}.sam > {0}.bam".format(prefix))
run_cmd("samtools sort {0}.bam > {0}.sorted.bam".format(prefix))
run_cmd("samtools index {0}.sorted.bam".format(prefix))