# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 14:27:47 2015

@author: stxjtd
"""

import multicell.simulation
import copy
import time
from multicell.simulation import print_flush
from openalea.container import property_topomesh
from vplants.meshing.property_topomesh_analysis import compute_topomesh_property
from openalea.container.utils                       import IdDict
import numpy as np

class SimulationPTM(multicell.simulation.Simulation):
    
    def __init__(self):
        multicell.simulation.Simulation.__init__(self)
        self.computed = set()
        
    def import_propertytopomesh(self, mesh):
        time_start = time.time()
        print_flush("Topomesh importation: started")
        self.set_mesh(mesh)
        pos = self.get_pos().values()
        self.set_pos(pos - np.mean(pos, axis=0))
        self.initialize_mesh_properties()
        print_flush("Topomesh importation: finished (%.2f s)" % (time.time() - time_start))
        
    def set_mesh(self, mesh):
        if type(mesh) is property_topomesh.PropertyTopomesh:
            self.mesh = mesh
        else:
            degree = 3
            
            # Manually build a PropertyTopomesh from a regular Topomesh
            # Dirty but necessary in the current version of openalea
            ptm = copy.deepcopy(mesh)
            ptm.__class__ = property_topomesh.PropertyTopomesh
            ptm._wisp_properties = [{} for d in xrange(degree+1)]
            ptm._interface_properties = [{} for d in xrange(degree+1)]
            ptm._topomesh_properties= {}
            ptm._interface = [None] + [IdDict(idgenerator="set") for i in xrange(degree)]

            self.mesh = ptm #property_topomesh.PropertyTopomesh(3, mesh)
    
    def get_pos(self):
        return self.mesh.wisp_property("barycenter", 0)
        
    def set_pos(self, pos):
        self.mesh.update_wisp_property("barycenter", 0, pos)
        for degree in xrange(1, 4):
            self.mesh._wisp_properties[degree].clear()
        self.computed = set()
    
    def compute_property(self, name, degree):
        compute_topomesh_property(self.mesh, name, degree)
        self.computed.add((name, degree))
        
    def get_property(self, name, degree, wid):
        return self.get_properties(name, degree)[wid]
    
    def get_properties(self, name, degree):
        if not (name, degree) in self.computed:
            self.compute_property(name, degree)
        return self.mesh.wisp_property(name, degree)
    
    def compute_surfaces(self):
        self.compute_property("area", 2)
        
    def get_surface(self, wid):
        return self.get_property("area", 2, wid)
    
    def get_surfaces(self):
        return self.get_properties("area", 2)
        
    def compute_barycenters(self):
        self.compute_property("barycenter", 3)
        
    def get_barycenter(self, wid):
        return self.get_property("barycenter", 3, wid)
    
    def get_barycenters(self):
        return self.get_properties("barycenter", 3)