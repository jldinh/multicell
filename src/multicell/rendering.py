# -*- coding: utf-8 -*-
"""
Created on Sun Jul 17 22:38:13 2016

@author: jl
"""

#from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import time
from multicell.utilities import print_flush
from openalea.container import topomesh_algo
import os

class MatplotlibRenderer():
    """A class handling the rendering of a Simulation object in PlantGL"""
    
    def __init__(self, sim, max_cmap=None, view_size=None, view=(None, None)):
        """
        Creates a PlantGLRenderer using a Simulation object
        
        Parameters
        ----------
            sim : Simulation
                The Simulation object containing the data to render.
        """
        
        self.sim = sim
        self.max_cmap = max_cmap
        self.view_size = view_size
        self.view = view

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
        
        outer_fids = [fid for fid in self.sim.mesh.wisps(2) if self.sim.mesh.nb_regions(2, fid) == 1]
        if name == None:
            face_values = [0. for fid in outer_fids]
            array = np.zeros(self.sim.n_cells)
            max_cmap = 1
        else:
            environment = self.sim.compute_environment()
            associated_cids = [self.sim.mesh.regions(2, fid).next() for fid in outer_fids]
            array = environment[name]
            face_values = [array.get_species(array.variables_list[0], cid) for cid in associated_cids]
            if self.max_cmap is not None:
                max_cmap = self.max_cmap
            elif max_percentile is None:
                max_cmap = np.max(array)
            else:
                max_cmap = np.percentile(array, max_percentile)
            if max_cmap == 0:
                max_cmap = 1
                
        sm = matplotlib.cm.ScalarMappable(matplotlib.colors.Normalize(vmin=0, vmax=max_cmap, clip=True), "jet")
        facecolors = [sm.to_rgba(x) for x in face_values]
        
        
        polys = [[self.sim.pos[pid] for pid in topomesh_algo.ordered_pids(self.sim.mesh, fid)] for fid in outer_fids]
        
        
        
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.set_aspect('equal')
        ax.view_init(*self.view)
        
        poly = Poly3DCollection(polys, facecolors=facecolors, linewidth=0.2)
        ax.add_collection3d(poly)
        
        min_coords = np.min(self.sim.pos.values(), axis=0)
        max_coords = np.max(self.sim.pos.values(), axis=0)
        if self.view_size is not None:
            max_half_amplitude = self.view_size / 2
        else:
            max_half_amplitude = max(max_coords - min_coords) / 2
        center = (max_coords + min_coords) / 2
        boundaries = np.tile(center[:, np.newaxis], (1, 2))
        boundaries[:, 0] -= max_half_amplitude
        boundaries[:, 1] += max_half_amplitude
        
        ax.set_xlabel('X')
        ax.set_xlim3d(boundaries[0][0], boundaries[0][1])
        ax.set_ylabel('Y')
        ax.set_ylim3d(boundaries[1][0], boundaries[1][1])
        ax.set_zlabel('Z')
        ax.set_zlim3d(boundaries[2][0], boundaries[2][1])
        
        plt.show()

        if save:
            directory = "figures"
            if not os.path.exists(directory):
                os.makedirs(directory)
            timestamp = time.strftime('%Y%m%d_%H%M%S') + "_" + str(time.time() * 1000 % 1)[2:]
            plt.savefig("figures/" + name + "-" + timestamp + ".png")
#            Viewer.frameGL.saveImage("figures/" + name + "-" + timestamp + ".png")
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
            environment = self.sim.compute_environment()
            for s in self.sim.intensive_cell_variables:
                concs = environment[s][0]
                print_flush("%s: from %s to %s" % (s, min(concs), max(concs)))
        if max_percentile <> None:
            print_flush("Max value displayed for %s: %s" % (name, max_cmap))
