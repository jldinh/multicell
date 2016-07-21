# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 15:22:21 2015

@author: jl
"""

import numpy as np
from multicell import utilities

def linear_growth(mesh, pos, coefficient):
    """Applies a homotety to a dictionary of coordinates.
    
    Parameters
    ----------
        mesh : Topomesh
            Not used in this algorithm
        pos : dict(int -> iterable)
            Dictionary (pid -> ndarray) of the tissue vertices
        coefficient : float or ndarray
            Scaling coefficient for the homothety
        
    Returns
    -------
        dict(int -> ndarray)
            dictionary (pid -> new position) of the vertices
    """
    
    utilities.check_pos(pos)
    scaling = np.array(coefficient)
    res = dict((pid, scaling * vec) for pid,vec in pos.iteritems())
    assert np.all(res.values() <> None)
    return res