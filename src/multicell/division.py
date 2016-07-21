# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 15:24:11 2015

@author: jl
"""

from openalea.tissueshape import cell_main_axes, cell_volume, centroid
from multicell import utilities
import numpy as np

def symmetrical_division(mesh, pos, cid):
    """Divide a cell, roughly in the middle.
    
    The division plane is perpendicular to the longest axis and goes through
    the barycenter of the cell.    
    """
    
    utilities.check_pos(pos)
    bary, V1, V2, V3 = cell_main_axes(mesh, pos, cid)
    return bary, V1

def division_with_normal(mesh, pos, cid, normal):
    bary = centroid(mesh, pos, 3, cid)
    return bary, np.array(normal)

def volume_trigger(mesh, pos, cid, volume_threshold):
    """Indicates if the volume of a cell is above the threshold.
    """
    
    utilities.check_pos(pos)
    volume = cell_volume(mesh, pos, cid)
    return volume > volume_threshold

def always_trigger(mesh, pos, cid):
    return True