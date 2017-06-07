__author__ = 'etseng@pacb.com'

from cupcake2.io.minimapIO import MiniReader
from cupcake2.ice2.preCluster import preClusterSet

NUM_SEED = 20000
NUM_INCOMING_SEQS = 1000000

def init_seed_S(minimap_filename, max_missed_5_len, max_missed_5_ratio, max_missed_3_len, max_missed_3_ratio):
    """
    seed minimap contains self-hits and redundant hits, so be sure to ignore them.
    """
    pCS = preClusterSet()

    for r in MiniReader(minimap_filename):
        if r.qID == r.sID or r.strand == '-': continue # ignore self and opp strand hits
        if r.qID > r.sID: continue # ignore redundant hits since this is all-for-all
        stat = r.characterize(max_missed_5_len, max_missed_5_ratio, max_missed_3_len, max_missed_3_ratio)
        # stat could be: match, q_contained, s_contained, partial
        # at this first stage, we only care about match
        if stat == 'match':
            pCS.add_seqid_match(r.qID, r.sID)
        elif stat == 'q_contained':
            pCS.add_tucked_match(r.qID, r.sID)
        elif stat == 's_contained':
            pCS.add_tucked_match(r.sID, r.qID)
        else:
            pass # do nothing if it's partial


