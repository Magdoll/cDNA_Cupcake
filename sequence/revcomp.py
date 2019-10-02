#!/usr/bin/env python

__version__ = '1.0'

import os, sys
from Bio import Seq

for seq in sys.argv[1:]:
    print(Seq.Seq(seq).reverse_complement())
