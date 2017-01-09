# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 16:05:13 2015

@author: jl
"""

import unittest
from multicell import simulation_builder, division
import numpy as np

class TestDivisionFunctions(unittest.TestCase):

    def test_symmetrical_division_returns_3D_vectors_in_3D_tissue(self):
        mesh, pos = simulation_builder.generate_cell_grid(1,1,1)
        cid = list(mesh.wisps(3))[0]
        point, normal = division.symmetrical_division(mesh, pos, cid)
        self.assertEqual(point.__class__, np.ndarray)
        self.assertEqual(normal.__class__, np.ndarray)
        self.assertEqual(len(point), 3)
        self.assertEqual(len(normal), 3)
    
    def test_volume_trigger_returns_True_if_volume_greater_than_threshold(self):
        mesh, pos = simulation_builder.generate_cell_grid(1,1,1)
        cid = list(mesh.wisps(3))[0]
        is_triggered = division.volume_trigger(mesh, pos, cid, 0.5)
        self.assertEqual(is_triggered, True)
        
    def test_volume_trigger_returns_False_if_volume_lower_than_threshold(self):
        mesh, pos = simulation_builder.generate_cell_grid(1,1,1)
        cid = list(mesh.wisps(3))[0]
        is_triggered = division.volume_trigger(mesh, pos, cid, 1.5)
        self.assertEqual(is_triggered, False)    
    
if __name__ == '__main__':
    unittest.main()