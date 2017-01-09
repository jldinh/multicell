# -*- coding: utf-8 -*-
"""
Created on Mon Dec 21 18:58:01 2015

@author: jean-louis
"""

import unittest
from multicell import simulation_builder, division, growth
import numpy as np

class TestSimulationMethods(unittest.TestCase):

    def test_cell_division_splits_a_cell_in_two_with_symmetrical_division_method(self):
        sim = simulation_builder.generate_cell_grid_sim(1,1,1)
        cid = list(sim.mesh.wisps(3))[0]
        point, normal = division.symmetrical_division(sim.mesh, sim.pos, cid)
        sim.cell_division(cid, point, normal)
        self.assertEqual(len(list(sim.mesh.wisps(3))), 2)
        
    def test_cell_division_splits_a_cell_in_two_with_parallel_plane(self):
        sim = simulation_builder.generate_cell_grid_sim(1,1,1)
        cid = list(sim.mesh.wisps(3))[0]
        point = np.array([0.5,0.5,0.5])
        normal = np.array([1,0,0])
        sim.cell_division(cid, point, normal)
        self.assertEqual(len(list(sim.mesh.wisps(3))), 2)
        
    def test_cell_division_only_attempts_to_split_relevant_edges_with_parallel_plane_and_no_noise(self):
        sim = simulation_builder.generate_cell_grid_sim(1,1,1,0.)
        cid = list(sim.mesh.wisps(3))[0]
        point = np.array([0.5,0.5,0.5])
        normal = np.array([1,0,0])
        sim.cell_division(cid, point, normal)
        self.assertEqual(len(list(sim.mesh.wisps(3))), 2)
    
    def test_minimal_simulation_works(self):
        sim = simulation_builder.generate_cell_grid_sim(1,1,1)
        sim.register_cell_variable("a")
        sim.initialize_cell_variables()
        sim.simulate()
    
    def test_division_works(self):
        sim = simulation_builder.generate_cell_grid_sim(1,1,1)
        sim.register_cell_variable("a")
        sim.initialize_cell_variables()
        sim.enable_division()
        sim.register_division_method(division.division_with_normal, {"normal": np.array([1,0,0])})
        sim.register_division_trigger(division.always_trigger)
        sim.simulate()
        self.assertEqual(sim.cids.size, 2)
        self.assertEqual(sim.n_cells, 2)
        self.assertEqual(sim.y.current().size, 2)
    
if __name__ == '__main__':
    unittest.main()