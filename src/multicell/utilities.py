# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 17:12:12 2015

@author: jl
"""

import sys
import numpy as np

def print_flush(msg):
    print msg
    sys.stdout.flush()

def check_pos(pos):
    if np.any([not isinstance(x, np.ndarray) for x in pos.values()]):
        raise Exception("Invalid pos: pos contains values that are not ndarrays.")