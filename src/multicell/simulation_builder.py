# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 11:34:49 2015

@author: jl
"""

import numpy as np
from openalea.tissueshape import tovec, grid_tissue
from multicell import simulation
        
def generate_cell_grid(x, y=0, z=0, noise_amplitude=1e-3):
    """Builds the representation of a grid tissue of the specified shape.
    
    Parameters
    ----------
        x : int
            Number of cells along the x axis.
        y : int
            Number of cells along the y axis, not used if 0. (default: 0)
        z : int
            Number of cells along the z axis, not used if 0. (default: 0)
        noise_amplitude : double
            Amplitude of the random noise applied to the position of vertices.
            The perturbations applied to each vertex are independent and
            follow a uniform distribution between -noise_amplitude/2. and
            noise_amplitude/2. (default: 1e-3)
    
    Returns
    -------
        Topomesh
            Topomesh of the grid tissue
        dict
            Dictionary (vertex ID -> position) of the tissue
    """
    
    shape = tuple([u for u in [x,y,z] if u >= 1])
    tissuedb = grid_tissue.regular_grid(shape)
    mesh = tissuedb.get_topology("mesh_id")
    pos = tovec(tissuedb.get_property("position") )
    bary = reduce(lambda x,y: x + y,pos.itervalues() ) / len(pos)
    pos = dict((pid,vec - bary + np.random.uniform(-noise_amplitude/2., noise_amplitude/2., len(shape))) for pid,vec in pos.iteritems())
    return mesh, pos
    
def generate_cell_grid_sim(x, y=0, z=0, noise_amplitude = 1e-3, sim_class=simulation.Simulation):
    """Builds a Simulation object with a grid-shaped tissue.
    
    Parameters
    ----------
        x : int
            Number of cells along the x axis
        y : int
            Number of cells along the y axis, not used if 0 (default: 0)
        z : int
            Number of cells along the z axis, not used if 0 (default: 0)
        noise_amplitude : double
            Amplitude of the random noise applied to the position of vertices.
            The perturbations applied to each vertex are independent and
            follow a uniform distribution between -noise_amplitude/2. and
            noise_amplitude/2. (default: 1e-3)
    
    Returns
    -------
        Simulation
            Simulation object containing the grid.
    """
    
    sim = sim_class()
    mesh, pos = generate_cell_grid(x, y, z, noise_amplitude)
    sim.import_topomesh(mesh, pos)
    return sim
