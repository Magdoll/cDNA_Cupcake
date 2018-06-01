import os, sys
import pysam

if not os.path.exists('polished.bam'):
    print >> sys.stderr, "polished.bam does not exist. Abort!"
    sys.exit(-1)

f = pysam.AlignmentFile('polished.bam', check_sq=False)
h = open('cluster_report.csv', 'w')
h.write('cluster_id,read_id,read_type\n')
for r in f:
    d = dict(r.tags)
    for x in d['im'].split(','):
        # massage ccs name into <movie>/<zmw>/<start_end_CCS>
        _movie, _zmw, _ccs = x.split('/')
        newid = "{0}/{1}/1_1000_CCS".format(_movie, _zmw)
        h.write("cb1_c{0},{1},FL\n".format(r.qname.split('/')[1], newid))
        
h.close()
print >> sys.stderr, "Output written to: cluster_report.csv"
