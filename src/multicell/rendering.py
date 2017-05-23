# -*- coding: utf-8 -*-
"""
Created on Sun Jul 17 22:38:13 2016

@author: jl
"""

#from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import time
from multicell.utilities import print_flush
from openalea.container import topomesh_algo
import os
from openalea.tissueshape import centroid

def unit_vector(vector):
    return vector / np.linalg.norm(vector)
    
def border_vector(corner, opposite, centroid):
#    adjacent = opposite - corner
    hypothenuse = centroid - corner
#    projection_length = np.dot(adjacent, hypothenuse) / np.linalg.norm(adjacent)
#    opposite_length = np.sqrt(np.linalg.norm(hypothenuse)**2 - projection_length**2)
#    result = hypothenuse / opposite_length
#    print result
#    return result
    return unit_vector(hypothenuse)

class MatplotlibRenderer(object):
    """A class handling the rendering of a Simulation object in Matplotlib"""
    
    def __init__(self, sim, max_cmap=None, view_size=None, view=(None, None), axes=True, clipping=((0,0,0), (0,0,0))):
        """
        Creates a Matplotlib renderer using a Simulation object
        
        Parameters
        ----------
            sim : Simulation
                The Simulation object containing the data to render.
        """
        
        self.sim = sim
        self.max_cmap = max_cmap
        self.view_size = view_size
        self.view = view
        self.axes = axes
        self.clipping = clipping

    def _render(self, name=None, save=False, max_percentile=None):
        """
        Low level method to render a tissue, colored by concentrations.
        
        Concentrations are taken from the table of concentrations of the 
        Simulation. Uses Jet as ColorMap.
        
        Parameters
        ----------
            name : string
                Name of the species whose concentrations must be
                rendered
            save : bool
                Whether to save a picture or not
        """
        if not any(self.clipping[1]):
            outer_fids = [fid for fid in self.sim.mesh.wisps(2) if self.sim.mesh.nb_regions(2, fid) == 1]
            non_clipped_cids = set(self.sim.mesh.wisps(3))
        else:
            non_clipped_cids = set()
            reference = np.array(self.clipping[0])
            normal = np.array(self.clipping[1])
            for cid in self.sim.mesh.wisps(3):
                vector = centroid(self.sim.mesh, self.sim.get_pos(), 3, cid) - reference
                if np.dot(vector, normal) > 0:
                    non_clipped_cids.add(cid)
            outer_fids = []
            for fid in self.sim.mesh.wisps(2):
                neighbors = set(self.sim.mesh.regions(2, fid))
                count = 0
                for cid in neighbors:
                    if cid in non_clipped_cids:
                        count += 1
                if count == 1:
                    outer_fids.append(fid)
        outer_eids = set()
        for fid in outer_fids:
            outer_eids.update(self.sim.mesh.borders(2, fid))
        cell_delimiting_eids = set()
        for eid in outer_eids:
            delimited_cells = set(self.sim.mesh.regions(1, eid, 2)) & non_clipped_cids
            if len(delimited_cells) > 1:
                cell_delimiting_eids.add(eid)
        if name == None:
            face_values = [0.5 for fid in outer_fids]
            array = np.zeros(self.sim.n_cells)
            max_cmap = 1
        else:
            environment = self.sim.compute_environment()
            associated_cids = [(set(self.sim.mesh.regions(2, fid)) & non_clipped_cids).pop() for fid in outer_fids]
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
        sm.set_array(face_values)
        facecolors = [sm.to_rgba(x) for x in face_values]
        
#        poly_pids = [[pid for pid in topomesh_algo.ordered_pids(self.sim.mesh, fid)] for fid in outer_fids]
        min_coords = np.min(self.sim.pos.values(), axis=0)
        max_coords = np.max(self.sim.pos.values(), axis=0)
        if self.view_size is not None:
            max_half_amplitude = self.view_size / 2
        else:
            max_half_amplitude = max(max_coords - min_coords) / 2
        shrinkage = 0.005 * max_half_amplitude
        cell_delimiting_pids = set()
        for eid in cell_delimiting_eids:
            cell_delimiting_pids.update(set(self.sim.mesh.borders(1, eid)))
        polys = []
        for fid in outer_fids:
            ordered_pids = list(topomesh_algo.ordered_pids(self.sim.mesh, fid))
            original_pos = [self.sim.pos[pid] for pid in ordered_pids]
            cid = (set(self.sim.mesh.regions(2, fid)) & non_clipped_cids).pop()
            bary = centroid(self.sim.mesh, self.sim.get_pos(), 3, cid)
            shrinked_pos = [original_pos[i] + (shrinkage * border_vector(original_pos[i], original_pos[(i+1)%len(ordered_pids)], bary)) for i in xrange(len(ordered_pids))]
            polys.append(shrinked_pos)
#            polys.append(original_pos)
        for eid in cell_delimiting_eids:
            delimited_cids = list(set(self.sim.mesh.regions(1, eid, 2)) & non_clipped_cids)
            centroids = [centroid(self.sim.mesh, self.sim.get_pos(), 3, cid) for cid in delimited_cids]
            edge_pids = list(self.sim.mesh.borders(1, eid))
            edge_poly = [self.sim.get_pos()[edge_pids[0]],
                         self.sim.get_pos()[edge_pids[0]] + shrinkage * border_vector(self.sim.get_pos()[edge_pids[0]], self.sim.get_pos()[edge_pids[1]], centroids[0]),
                         self.sim.get_pos()[edge_pids[1]] + shrinkage * border_vector(self.sim.get_pos()[edge_pids[1]], self.sim.get_pos()[edge_pids[0]], centroids[0]),
                         self.sim.get_pos()[edge_pids[1]],
                         self.sim.get_pos()[edge_pids[1]] + shrinkage * border_vector(self.sim.get_pos()[edge_pids[1]], self.sim.get_pos()[edge_pids[0]], centroids[1]),
                         self.sim.get_pos()[edge_pids[0]] + shrinkage * border_vector(self.sim.get_pos()[edge_pids[0]], self.sim.get_pos()[edge_pids[1]], centroids[1])]
            polys.append(edge_poly)
            facecolors.append((0,0,0,1))
#            
#        cell_delimiting_eids_as_pid_sets = set()
#        for eid in cell_delimiting_eids:
#            edge = tuple(self.sim.mesh.borders(1, eid))
#            cell_delimiting_eids_as_pid_sets.update([edge, tuple(reversed(edge))])
#        linewidths = []
#        for poly in poly_pids:
#            if np.any(poly[-1] != poly[0]):
#                poly.append(poly[0])
#            n = len(poly)
#            for i in xrange(n-1):
#                a = poly[i]
#                if i < n - 1:
#                    b = poly[i+1]
#                else:
#                    b = poly[0]
#                edge = (a, b)
##                print edge
#                if edge in cell_delimiting_eids_as_pid_sets:
#                    linewidths.append(1)
#                else:
#                    linewidths.append(0)
#        segments = [[self.sim.pos[pid] for pid in self.sim.mesh.borders(1, eid)] for eid in cell_delimiting_eids]
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.set_aspect('equal')
        ax.view_init(*self.view)
        if not self.axes:
            ax.set_axis_off()
        
#        print linewidths
#        print np.mean(linewidths)
#        linewidths+=[0.5]*10
        poly = Poly3DCollection(polys, facecolors=facecolors, linewidth=0)
#        edge_collection = Line3DCollection(segments, colors=((0,0,0,1),), linewidth=1.)
        ax.add_collection3d(poly)
#        ax.add_collection3d(edge_collection)
        
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
        
        if name is not None:
            fig.colorbar(sm, shrink=0.5, aspect=10)
            ax.set_title(name)
        
#        plt.subplots_adjust(left=0., right=1., top=1., bottom=0.)
        plt.tight_layout()
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
        print_flush("Time point: %s" % self.sim.current_t)       
        if name <> None:
            environment = self.sim.compute_environment()
            for s in self.sim.intensive_cell_variables:
                concs = environment[s][0]
                print_flush("%s: from %s to %s" % (s, min(concs), max(concs)))
        if max_percentile <> None:
            print_flush("Max value displayed for %s: %s" % (name, max_cmap))
