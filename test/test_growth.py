# -*- coding: utf-8 -*-
"""
Created on Mon Dec 21 18:11:58 2015

@author: jean-louis
"""

import unittest
from multicell import simulation_builder, growth
import numpy as np

class TestGrowthFunctions(unittest.TestCase):

    def test_growth_returns_dictionary_of_ints_and_ndarrays(self):
        mesh, pos = simulation_builder.generate_cell_grid(1,1,1)
        pos2 = growth.linear_growth(mesh, pos, 1.1)
        for cid, coords in pos2.items():
            self.assertIsInstance(cid, int)
            self.assertIsInstance(coords, np.ndarray)
            self.assertEqual(len(coords), 3)
        
if __name__ == '__main__':
    unittest.main()