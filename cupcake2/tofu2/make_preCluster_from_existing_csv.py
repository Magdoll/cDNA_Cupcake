
"""
If a sequence output CSV is already available (via some other preClustering method, etc),
generate the preCluster_out/ directory and other files so ICE can be run.

seqid,stat
m54006_170206_215027/39191167/4082_54_CCS,11800
m54086_170204_081430/54592493/1354_58_CCS,5900

seqid --- input FLNC read id
stat --- can be (cid) or "orphan"
"""
import os, sys
from csv import DictReader

from cupcake.io.SeqReaders import LazyFastaReader
from cupcake2.ice2.preCluster import preClusterSet2, preCluster
from cupcake2.tofu2.run_preCluster import detect_PCR_chimeras
import cupcake2.io.FileIO as FileIO

def read_seq_csv(csv_filename):
    # sanity check that "seqid" and "stat" are two valid column headers
    header_checked = False
    orphans = set()
    pCS = preClusterSet2()

    for r in DictReader(open(csv_filename), delimiter=','):
        if not header_checked:
            if 'seqid' not in r or 'stat' not in r:
                print >> sys.stderr, "{0} must have the fields 'seqid' and 'stat'! Abort".format(csv_filename)
                sys.exit(-1)
        header_checked = True
        if r['stat']=='orphan':
            orphans.add(r['seqid'])
        else:
            cid = int(r['stat'])
            if cid not in pCS.S: pCS.S[cid] = preCluster(cid=cid)
            pCS.add_seqid_to_cluster_by_cid(r['seqid'], cid)
    return pCS, orphans


def main(fasta_filename, csv_filename):
    d = LazyFastaReader(fasta_filename)
    pCS, orphans = read_seq_csv(csv_filename)

    # detect PCR chimeras from orphans
    chimeras = detect_PCR_chimeras(orphans, d)
    orphans = orphans.difference(chimeras)

    FileIO.write_seqids_to_fasta(orphans, "preCluster_out.orphans.fasta", d)
    FileIO.write_seqids_to_fasta(chimeras, "preCluster_out.chimeras.fasta", d)


    infof = open('preCluster.cluster_info.csv', 'w')
    infof.write("cluster,size\n")
    # write out a directory per preCluster cid in preCluster_out/<cid>
    # Liz note: right now, write out even directories with just 1 sequence
    # (we know they have "tucked" support, so can run Partial/Arrow on it)
    #singlef = open("preCluster_out.singles.fasta", 'w')
    for cid in pCS.S:
    #    if pCS.S[cid].size == 1:
    #        r = d[pCS.S[cid].members[0]]
    #        singlef.write(">{0}\n{1}\n".format(r.id, r.seq))
    #    else:
        #print >> sys.stderr, "cid", cid
        if True:
            dirname = os.path.join("preCluster_out", str(cid))
            os.makedirs(dirname)
            file = os.path.join(dirname, 'isoseq_flnc.fasta')
            FileIO.write_seqids_to_fasta(pCS.S[cid].members, file, d)
        infof.write("{0},{1}\n".format(cid, len(pCS.S[cid].members)))
        #print cid, len(pCS.S[cid].members)
    #singlef.close()
    infof.close()

if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser("Make preCluster bins based on external cluster csv result")
    parser.add_argument("fasta_filename")
    parser.add_argument("cluster_csv")

    args = parser.parse_args()
    main(args.fasta_filename, args.cluster_csv)