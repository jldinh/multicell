# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 16:18:16 2015

@author: jl
"""

import unittest
from multicell import simulation_builder, simulation

class TestSimulationBuilderMethods(unittest.TestCase):

    def test_generate_1D_tissue(self):
        mesh, _ = simulation_builder.generate_cell_grid(10)
        self.assertEqual(len(list(mesh.wisps(1))), 10)

    def test_generate_2D_tissue(self):
        mesh, _ = simulation_builder.generate_cell_grid(10, 10)
        self.assertEqual(len(list(mesh.wisps(2))), 100)

    def test_generate_3D_tissue(self):
        mesh, _ = simulation_builder.generate_cell_grid(10, 10, 10)
        self.assertEqual(len(list(mesh.wisps(3))), 1000)
    
    def test_generate_cell_grid_sim_returns_Simulation_object(self):
        sim = simulation_builder.generate_cell_grid_sim(1,1,1)
        self.assertEqual(sim.__class__, simulation.Simulation)
        
if __name__ == '__main__':
    unittest.main()