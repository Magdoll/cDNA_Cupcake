__author__ = 'etseng@pacb.com'

import os, sys, subprocess

"""
/usr/bin/time -v minimap2 -f 0.00001 -L 100 -t 30 --eqx --MD seed0.fasta seed0.fasta > seed0.fasta.f00001.minimap 2> seed0.fasta.f00001.minimap.log
"""

CMD_MINIMAP = \
"/usr/bin/time -v minimap2 {a} -H -f 0.00001 -N 2000000 -m 20 -X -c -t {c} {t} {q} > {o} 2> {o}.log"
OUTPUT_MINIMAP = "{q}.{t}.f00001.minimap"

def run_minimap(query_filename, target_filename, cpus, sam_output=False):
    output_filename = OUTPUT_MINIMAP.format(q=query_filename, t=target_filename)
    cmd = CMD_MINIMAP.format(q=query_filename, t=target_filename, \
                             o=output_filename, c=cpus, \
                             a="-a" if sam_output else "")
    run_cmd(cmd)
    return output_filename

def run_cmd(cmd):
    if subprocess.check_call(cmd, shell=True)!=0:
        print("ERROR RUNNING CMD:", cmd, file=sys.stderr)
        sys.exit(-1)


