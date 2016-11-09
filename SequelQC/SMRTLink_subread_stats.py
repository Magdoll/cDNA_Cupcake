__author__ = 'etseng@pacb.com'


"""
Getting # of ZMWs and # of subreads and subread lengths from subread XML file.
Requires pbcore from SA3.x! So putting this in a separate script for dev.

Note: also "approximates" ZMW length by using last subread end
"""

import os, sys
from collections import defaultdict
from pbcore.io import SubreadSet

def get_subread_ZMW_stats(subread_xml, report):
    """
    Fills a dict with:
    'numZMW' --- number of sequencing ZMWs
    'numSubread' -- number of subreads
    'avgZMWlen' -- approximated average ZMW length
    'avgSubreadlen' --- average subread length
    """
    subread_lens = []
    zmw_lens = defaultdict(lambda: 0)

    ds = SubreadSet(subread_xml)
    for rr in ds.resourceReaders():
        for zmw, qStart, qEnd in zip(rr.holeNumber, rr.qStart, rr.qEnd):
            subread_lens.append(qEnd-qStart)
            zmw_lens[zmw] = max(zmw_lens[zmw], qEnd)

    report['numZMW'] = len(zmw_lens)
    report['numSubread'] = len(subread_lens)
    report['avgZMWlen'] = int(sum(zmw_lens.itervalues())*1./len(zmw_lens))
    report['avgSubreadlen'] = int(sum(subread_lens)*1./len(subread_lens))

    #return zmw_lens

