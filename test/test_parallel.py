## -*- coding: utf-8 -*-
#"""
#Created on Mon Dec 21 18:58:01 2015
#
#@author: jean-louis
#"""
#
#import unittest
#from multicell import simulation_builder, parallel
#import numpy as np
#import multiprocessing.sharedctypes
#import ctypes
#
#class TestParallel(unittest.TestCase):
#    
#    def test_rawarray_attribute_of_concentration_table_is_original_rawarray(self):
#        cids = [0,1]
#        rawarray = multiprocessing.sharedctypes.RawArray(ctypes.c_double, 2)
#        ct = parallel.ConcentrationTableMultiprocessing(["A"], cids, cids, rawarray)
#        self.assertIs(rawarray, ct.rawarray)
#    
#    def test_concentration_table_built_from_other_table_s_rawarray_has_the_same_rawarray_attribute(self):
#        cids = [0,1,2,3]
#        ct1 = parallel.ConcentrationTableMultiprocessing(["A"], cids)
#        ct2 = parallel.ConcentrationTableMultiprocessing(["A"], cids, cids, ct1.rawarray)
#        self.assertIs(ct1.rawarray, ct2.rawarray)
#    
#    def test_modifying_rawarray_from_standard_array_affects_concentration_table(self):
#        cids = [0,1]
#        rawarray = multiprocessing.sharedctypes.RawArray(ctypes.c_double, 2)
#        ct = parallel.ConcentrationTableMultiprocessing(["A"], cids, cids, rawarray)
#        a = np.frombuffer(ct.rawarray)
#        a[...] = [3., 14.]
#        print ct
#        self.assertTrue(np.all(ct[0] == a))
#    
#    def test_modifying_concentration_table_writes_to_the_original_rawarray(self):
#        cids = [0,1]
#        rawarray = multiprocessing.sharedctypes.RawArray(ctypes.c_double, 2)
#        ct = parallel.ConcentrationTableMultiprocessing(["A"], cids, cids, rawarray)
#        values = np.array([3., 14.])
#        ct.set_species("A", values)
#        a = np.frombuffer(ct.rawarray)
#        print a
#        self.assertTrue(np.all(ct.get_species("A") == a))
#    
#    def test_modifying_concentration_table_writes_to_its_rawarray(self):
#        cids = [0,1]
#        ct = parallel.ConcentrationTableMultiprocessing(["A"], cids)
#        values = np.array([3., 14.])
#        ct.set_species("A", values)
#        a = np.frombuffer(ct.rawarray)
#        print a
#        self.assertTrue(np.all(ct.get_species("A") == a))
#        
#    def test_modifying_rawarray_from_two_concentration_tables_affect_both_tables(self):
#        cids = [0,1]
#        rawarray = multiprocessing.sharedctypes.RawArray(ctypes.c_double, 2)
#        ct1 = parallel.ConcentrationTableMultiprocessing(["A"], cids, cids, rawarray)
#        values = np.array([3., 14.])
#        ct1.set_species("A", values)
#        ct2 = parallel.ConcentrationTableMultiprocessing(["A"], cids, cids, rawarray)
#        print ct2
#        self.assertTrue(np.all(ct2.get_species("A") == values))
#    
#    def test_concentration_table_can_be_written_from_different_objects(self):
#        cids = [0,1,2,3]
#        ct = parallel.ConcentrationTableMultiprocessing(["A"], cids)
#        self.assertEqual(ct.shape, (1,len(cids)))
#        ct1 = parallel.ConcentrationTableMultiprocessing(["A"], [1], cids, ct.rawarray)
#        self.assertEqual(ct1.shape, (1,1))
#        ct2 = parallel.ConcentrationTableMultiprocessing(["A"], [2], cids, ct.rawarray)
#        a3 = np.frombuffer(ct.rawarray)
#        a3[3] = 3
#        ct1.set_species("A", 1)
#        ct2.set_species("A", np.array([2]))
#        ct_bis = parallel.ConcentrationTableMultiprocessing(["A"], cids, cids, ct.rawarray)
#        print ct
#        print ct1
#        print ct2
#        print ct_bis
#        self.assertEqual(ct.get_species("A", 0), 0)
#        self.assertEqual(ct.get_species("A", 1), 1)
#        self.assertEqual(ct.get_species("A", 2), 2)
#        self.assertEqual(ct.get_species("A", 3), 3)
#
##    def test_cids_are_partitioned_correctly(self):
##        sim = simulation_builder.generate_cell_grid_sim(1, 1, 50, sim_class=parallel.SimulationMultiprocessing)
##        sim.c
#
#    def test_derivative_is_not_all_zeros(self):
#        sim = simulation_builder.generate_cell_grid_sim(50, 1, 1, sim_class=parallel.SimulationMultiprocessing)
#        sim.register_cell_variable("A")
#        sim.set_ODE("A", "1")
#        sim.initialize_cell_variables()
#        sim.compile_ODEs()
#        d = sim.derivative(sim.y.current().as_1d_array(), 0)
#        self.assertFalse(np.all(d == 0))
#        self.assertTrue(np.all(d == 1))
#        
#    def test_constant_arrays_are_split_correctly(self):
#        sim = simulation_builder.generate_cell_grid_sim(50, 1, 1, sim_class=parallel.SimulationMultiprocessing)
#        sim.register_cell_variable("A")
#        sim.register_cell_variable("B")
#        sim.set_constant("k", np.ones(50), (True,))
#        sim.register_computed_variable("k1", "1*k")
#        sim.set_ODE("A", "k")
#        sim.set_ODE("B", "k1")
#        sim.initialize_cell_variables()
#        sim.compile_ODEs()
#        d = sim.derivative(sim.y.current().as_1d_array(), 0)
#        self.assertFalse(np.all(d == 0))
#        self.assertTrue(np.all(d == 1))
#    
#if __name__ == '__main__':
#    unittest.main()