#!/usr/bin/env python
__author__ = 'etseng@pacb.com'

import os, sys, glob, datetime
import numpy as np
import json
from Bio import SeqIO
import SMRTLink_subread_stats as subread_stats

"""
Collect stats from a completed SMRTLink Iso-Seq jobs


CCS Report: pbreports.tasks.ccs_report-0/ccs_report.json
"""

LOCATION_ccs_report = 'tasks/pbreports.tasks.ccs_report-0/ccs_report.json'
LOCATION_flnc_report = 'tasks/pbreports.tasks.isoseq_classify-0/isoseq_classify_report.json'
LOCATION_flnc_pattern = 'tasks/pbtranscript.tasks.classify-*/isoseq_flnc.fasta'
LOCATION_cluster_report = 'tasks/pbreports.tasks.isoseq_cluster-0/isoseq_cluster_report.json'

LOCATION_smrtpipe_log = 'logs/pbsmrtpipe.log'

#LOCATION_ccs_stdout = 'tasks/pbreports.tasks.ccs_report-0/stdout'
#LOCATION_classify_stdout = 'tasks/pbreports.tasks.isoseq_classify-0/stdout'
#LOCATION_sep_flnc_stdout = 'tasks/pbtranscript.tasks.separate_flnc-0/stdout'

def get_subread_xml(smrtlink_dir):
    d = eval(os.popen("grep \"'eid_subread'\" {0}/job.stdout".format(smrtlink_dir)).read().strip())
    return d['eid_subread']

def read_ccs_report(smrtlink_dir, report):
    file = os.path.join(smrtlink_dir, LOCATION_ccs_report)
    d = json.load(open(file))

    for x in d['attributes']:
        if x['id'] == 'ccs.number_of_ccs_reads':
            report['numCCSread'] = int(x['value'])
        elif x['id'] == 'ccs.mean_ccs_readlength':
            report['avgCCSlen'] = int(x['value'])


def read_flnc_report(smrtlink_dir, report):
    file = os.path.join(smrtlink_dir, LOCATION_flnc_report)
    d = json.load(open(file))

    for x in d['attributes']:
        if x['id'] == 'isoseq_classify.num_nfl':
            report['numNFL'] = int(x['value'])
        elif x['id'] == 'isoseq_classify.num_flnc':
            report['numFLNC'] = int(x['value'])
        elif x['id'] == 'isoseq_classify.avg_flnc_len':
            report['avgFLNClen'] = int(x['value'])


def collect_flnc(smrtlink_dir, outdir):
    """
    Collect the chunked isoseq_flnc.fasta into one also collect lengths.
    Return: list of FLNC lengths
    """
    f = open(os.path.join(outdir, 'isoseq_flnc.fasta'), 'w')
    flnc_lens = []
    for file in glob.glob(os.path.join(smrtlink_dir, LOCATION_flnc_pattern)):
        for r in SeqIO.parse(open(file), 'fasta'):
            f.write(">{0}\n{1}\n".format(r.id, r.seq))
            flnc_lens.append(len(r.seq))
    f.close()
    return flnc_lens


def read_cluster_report(smrtlink_dir, report):
    file = os.path.join(smrtlink_dir, LOCATION_cluster_report)
    d = json.load(open(file))

    for x in d['attributes']:
        if x['id'] == 'isoseq_cluster.num_polished_hq_isoforms':
            report['numHQ'] = int(x['value'])


def collect_runtimes(smrtlink_dir, report):
    """
    Collect runtimes for: CCS, classify, cluster

    Don't precisely know the end points of CCS and classify -- instead approximate.
    CCS -- start of CCS to beginning of classify
    classify --- start of classify to beginning of cluster
    cluster --- start of cluster to end of task (including writing reports and IO)
    """
    ccs_start_t = None
    classify_start_t = None
    cluster_start_t = None
    partial_start_t = None
    polish_start_t = None
    task_end_t = None
    # collect CCS runtime
    for line in open(os.path.join(smrtlink_dir, LOCATION_smrtpipe_log)):
        # for CCS -- record the first CCS worker submitted time
        if line.startswith('[INFO]') and line.find('Starting worker pbccs.tasks.ccs') > 0 and ccs_start_t is None:
            raw = line.strip().split()
            ccs_start_t = datetime.datetime.strptime(raw[1]+' '+raw[2].split(',')[0], '%Y-%m-%d %H:%M:%S')
        elif line.startswith('[INFO]') and line.find('Starting worker pbtranscript.tasks.classify') > 0 and classify_start_t is None:
            raw = line.strip().split()
            classify_start_t = datetime.datetime.strptime(raw[1]+' '+raw[2].split(',')[0], '%Y-%m-%d %H:%M:%S')
        elif line.startswith('[INFO]') and line.find('Starting worker pbtranscript.tasks.cluster_bins') > 0 and cluster_start_t is None:
            raw = line.strip().split()
            cluster_start_t = datetime.datetime.strptime(raw[1]+' '+raw[2].split(',')[0], '%Y-%m-%d %H:%M:%S')
        elif line.startswith('[INFO]') and line.find('Starting worker pbtranscript.tasks.ice_partial_cluster_bins') > 0 and partial_start_t is None:
            raw = line.strip().split()
            partial_start_t = datetime.datetime.strptime(raw[1]+' '+raw[2].split(',')[0], '%Y-%m-%d %H:%M:%S')
        elif line.startswith('[INFO]') and line.find('Starting worker pbtranscript.tasks.ice_polish_cluster_bins') > 0 and polish_start_t is None:
            raw = line.strip().split()
            polish_start_t = datetime.datetime.strptime(raw[1]+' '+raw[2].split(',')[0], '%Y-%m-%d %H:%M:%S')
        elif line.startswith('[INFO]') and line.find('Starting worker pbtranscript.tasks.combine_cluster_bins') > 0:
            raw = line.strip().split()
            polish_end_t = datetime.datetime.strptime(raw[1]+' '+raw[2].split(',')[0], '%Y-%m-%d %H:%M:%S')
        elif line.startswith('[INFO]') and line.find('Completed execution pbsmrtpipe') > 0:
            raw = line.strip().split()
            task_end_t = datetime.datetime.strptime(raw[1]+' '+raw[2].split(',')[0], '%Y-%m-%d %H:%M:%S')

    report['runtimeCCS'] = (classify_start_t - ccs_start_t).total_seconds()
    report['runtimeFLNC'] = (cluster_start_t - classify_start_t).total_seconds()
    report['runtimeICE'] = (partial_start_t - cluster_start_t).total_seconds()
    report['runtimePolish'] = (polish_end_t - partial_start_t).total_seconds()
    report['runtimeTotal'] = (task_end_t - ccs_start_t).total_seconds()


def collect_RCvalidate_result(report):
    report['numIso'] = 'NA'

    iso_lens = []

    if os.path.exists('tmp'): # evidence that RCvalidate was already run
        for r in SeqIO.parse(open('tmp/touse.rep.fq'), 'fastq'):
            iso_lens.append(len(r.seq))
        report['numIso'] = len(iso_lens)

    return iso_lens


def gather_all_reports(smrtlink_dir, sample, out_dir='.'):
    report = {}

    subread_xml = get_subread_xml(smrtlink_dir)
    subread_stats.get_subread_ZMW_stats(subread_xml, report)

    read_ccs_report(smrtlink_dir, report)
    read_flnc_report(smrtlink_dir, report)
    read_cluster_report(smrtlink_dir, report)

    flnc_lens = collect_flnc(smrtlink_dir, out_dir)

    collect_runtimes(smrtlink_dir, report)

    iso_lens = collect_RCvalidate_result(report)


    with open('smrtlink_isoseq.report.txt', 'w') as f:
        f.write("Sample,Cells,ZMWs,Avg. ZMW length,Subreads,Avg. Subread length,CCS reads,Avg. CCS length,FLNC reads,FLNC %,Avg. FLNC length,FLNC Length Range,CCS runtime,Classify runtime,ICE runtime,Polish runtime,Total runtime")
        f.write(",HQ isoforms,Unique isoforms,Isoform Lengths\n")
        f.write("{sample},{cells},{zmw},{zmwlen},{sub},{sublen},{CCS},{CCSlen},{FLNC},{FLperc},{FLNClen},{FLNCrange},{runCCS},{runFLNC},{runICE},{runPolish},{runTotal}".format(\
            sample=sample, cells=1, zmw=report['numZMW'], zmwlen=report['avgZMWlen'],\
            sub=report['numSubread'], sublen=report['avgSubreadlen'],\
            CCS=report['numCCSread'], CCSlen=report['avgCCSlen'], FLNC=report['numFLNC'],\
            FLperc=report['numFLNC']*1./report['numCCSread'], FLNClen=report['avgFLNClen'],\
            FLNCrange="{0:.0f}-{1:.0f}".format(np.percentile(flnc_lens, 25), np.percentile(flnc_lens, 75)),\
            runCCS=report['runtimeCCS'], runFLNC=report['runtimeFLNC'],\
            runICE=report['runtimeICE'], runPolish=report['runtimePolish'],\
            runTotal=report['runtimeTotal']))
        f.write(",{hq},{iso},{isorange}\n".format(hq=report['numHQ'], \
            iso=report['numIso'], \
            isorange="{0:.0f}-{1:.0f}".format(np.percentile(iso_lens, 25), np.percentile(iso_lens, 75))))



if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser("Get stats from SMRTLink Iso-Seq completed job")
    parser.add_argument("smrtlink_dir")

    args = parser.parse_args()

    sample = os.path.basename(os.getcwd())
    gather_all_reports(args.smrtlink_dir, sample)
