# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 11:21:08 2015

@author: jl
"""

import numpy as np
import time
from openalea.plantgl.all import Scene, Viewer
from openalea.plantgl.ext.color import JetMap
import guillaume.topomesh_display
from multicell.utilities import print_flush

class PlantGLRenderer():
    """A class handling the rendering of a Simulation object in PlantGL"""
    
    def __init__(self, sim):
        """
        Creates a PlantGLRenderer using a Simulation object
        
        Parameters
        ----------
            sim : Simulation
                The Simulation object containing the data to render.
        """
        
        self.sim = sim

    def _render(self, name=None, save=False, max_percentile=None):
        """
        Low level method to render a tissue, colored by concentrations.
        
        Concentrations are taken from the table of concentrations of the 
        Simulation. Uses JetMap as ColorMap.
        
        Parameters
        ----------
            name : string
                Name of the species whose concentrations must be
                rendered
            save : bool
                Whether to save a picture or not
        """
        
        if name == None:
            array = np.zeros(self.sim.n_cells)
            max_cmap = 0
        else:
            array = self.sim.y.get_species(name) / self.sim.dilution_volumes.as_1d_array()
            if max_percentile is None:
                max_cmap = np.max(array)
            else:
                max_cmap = np.percentile(array, max_percentile)
                array = np.where(array <= max_cmap, array, max_cmap * np.ones(array.shape))
        prop = dict(zip(self.sim.cids, array))
        scene = Scene()
        scene += (guillaume
                    .topomesh_display
                    .drawPropertyTopomesh(self.sim.mesh, 
                                          self.sim.get_pos(),
                                          3,
                                          prop,
                                          JetMap,
                                          color_range=[0., max(max_cmap, 1e-6)],
                                          coef=0.99,
                                          transparent_min=False))
        Viewer.display(scene)
        if save:
            timestamp = time.strftime('%Y%m%d_%H%M%S') + "_" + str(time.time() * 1000 % 1)[2:]
            Viewer.frameGL.saveImage("figures/" + name + "-" + timestamp + ".png")
        return max_cmap

    def display(self, name=None, save=False, max_percentile=None):
        """
        Display the tissue with cells colored based on their concentrations.
        
        Also displays the lowest and highest concentrations in the tissue for
        all species
        
        Parameters
        ----------
        name : string
            Name of the species to display
        save : bool
            Whether a picture should be saved or not
        """
        
        max_cmap = self._render(name, save=save, max_percentile=max_percentile)
        if name <> None:
            for s in self.sim.names_species:
                concs = self.sim.y.get_species(s) / self.sim.dilution_volumes.as_1d_array()
                print_flush("%s: from %s to %s" % (s, min(concs), max(concs)))
        print_flush("Max value displayed for %s: %s" % (name, max_cmap))