# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 18:32:43 2016

@author: jl
"""

import ctypes
import multicell
import multiprocessing
import multiprocessing.sharedctypes
import numpy as np
import math
import itertools
import multiprocess
import traceback
import sys

class ConcentrationTableMultiprocessing(multicell.concentration_table.Concentration_Table):
    
    def __new__(cls, variables_list, cids, all_cids=None, rawarray=None):
        if rawarray is None:
#            print "rawarray is None"
            shared_mem = multiprocessing.sharedctypes.RawArray(ctypes.c_double, len(variables_list) * len(cids))
        else:
            shared_mem = rawarray
            
        if all_cids is None:
            acids = cids
        else:
            acids = all_cids
            
        array0 = np.frombuffer(shared_mem, dtype=ctypes.c_double)
        
        if rawarray is None:
            array0.fill(0)
        
#        print array.shape, (len(variables_list), len(acids))
        
        array0.shape = (len(variables_list), len(acids))
        
        slice_ = subset_slice(cids, acids)
        
        array = array0[:, slice_].view(cls)
        
        assert np.may_share_memory(array, array0)
        
        array.variables_list = variables_list
        array.variables_dict = dict((name, i) for (i, name) in enumerate(variables_list))
        array.cids_array = np.array(cids)
#        array.lcids_array = np.array(acids)
        array.dict_cids = dict((cid, i) for (i, cid) in enumerate(cids))
        array.n_cells = len(cids)
        array.rawarray = shared_mem
#        array.full_array = array0
        return array
    
    def __array_finalize__(self, a):
        super(ConcentrationTableMultiprocessing, self).__array_finalize__(a)
        if a is None:
            return
#        self.lcids_array = getattr(a, 'lcids_array', None)
        self.rawarray = getattr(a, 'rawarray', None)
#        self.full_array = getattr(a, 'full_array', None)
        
class SimulationMultiprocessing(multicell.simulation_ptm.SimulationPTM):
    
    def __init__(self, *args, **kwargs):
        super(SimulationMultiprocessing, self).__init__(*args, **kwargs)
        self.variable_split_axes = {}
    
    def initialize_mesh_properties(self):
        super(SimulationMultiprocessing, self).initialize_mesh_properties()
        self.sort_cids()
    
    def sort_cids(self, axis=0, n_bins=4):
        cids = self.get_barycenters().keys()
        barycenters = self.get_barycenters().values()
        percentile_cut_offs = np.linspace(0, 100, n_bins + 1)
        cut_offs = [np.percentile(barycenters[:, axis], x) for x in percentile_cut_offs]
        cut_offs[-1] += 1
        associated_bins = np.digitize(barycenters[:, axis], cut_offs) - 1
#        print zip(cids, associated_bins)
        
        bins = []
        n_sub_bins = 3 * n_bins - 2
        for i in xrange(n_sub_bins):
            bins.append(set())
            if i % 3 == 0:
                bins[i].update(cids[associated_bins == i / 3])
        
        try:        
            for i in xrange(0, n_sub_bins - 3, 3):
                for cid in bins[i]:
                    for ocid in self.mesh.border_neighbors(3, cid):
                        if ocid in bins[i + 3]:
#                            print i, bins[i], cid, ocid
                            bins[i + 1].add(cid)
                            bins[i + 2].add(ocid)
        except KeyError:
            raise KeyError("Slicing too thin. Decrease the number of bins or change slicing axis.")

        for i in xrange(1, n_sub_bins - 2, 3):
            bins[i - 1].difference_update(bins[i])
            bins[i + 2].difference_update(bins[i + 1])
                    
        
        associated_sub_bins = np.array(list(itertools.chain(*[[i] * len(b) for i,b in enumerate(bins)])))
        
        order = np.argsort(associated_sub_bins)
        print len(associated_sub_bins), len(cids)
        assert len(associated_sub_bins) == len(cids)
        self.cids = cids[order]
        self.bins = associated_sub_bins[order]
        self.n_bins = n_bins
        
    def compile_ODEs(self):
        """
        Compiles the expressions of ODEs and prepares a function to feed into
        the solver.
        """
        global dydt, ct, pool
        
        if hasattr(self, "pool"):
            self.pool.close()
        
        self.compute_dependencies()
        self.compute_Jacobian()
        self.derivative_components = {}
        for name in self.names_species:
            self.derivative_components[name] = compile(self.ODEs[name], "dydt_" + name, "eval")
        
#        n_processes = 8
#        barycenters = self.get_barycenters()
#        inner_cids = []
#        cids = []
#        for i in xrange(n_processes):
#            cids.append(set())
#            inner_cids.append(set())
#            
#        for cid, coords in barycenters.items():
#            bin_ = np.sum((coords > np.median(barycenters.values(), axis=0)) * (2 ** np.arange(3)))
#            inner_cids[bin_].add(cid)
#            cids[bin_].add(cid)
#            cids[bin_].update(self.mesh.border_neighbors(3, cid))
        
        n_sub_bins = 3 * self.n_bins - 2
        cids = [self.cids[np.logical_and(i - 2 <= self.bins, self.bins<= i + 2)] for i in xrange(0, n_sub_bins, 3)]
        inner_cids = [self.cids[np.logical_and(i - 1 <= self.bins, self.bins<= i + 1)] for i in xrange(0, n_sub_bins, 3)]
#        print self.bins
#        print cids
#        print inner_cids
        
        dydt = multicell.parallel.ConcentrationTableMultiprocessing(self.names_species, self.cids)   
        ct = multicell.parallel.ConcentrationTableMultiprocessing(self.names_species, self.cids)
        
        pool = multiprocess.Pool(initializer=init, initargs=(dydt.rawarray, ct.rawarray, self))        
        self.pool = pool
        
        def derivative(y, t):
            global dydt, ct, pool
            # Initialization of the derivative vector
            dydt.fill(0)
            ct.import_values(y)
            ct *= (ct>0)
            
            # multiprocessing
            
            pool.map(work, [(t, cids[i], inner_cids[i]) for i in xrange(self.n_bins)])
#            print dydt
#            pool.join()
            
            result = dydt.as_1d_array()

            # Test
            #print len(result), len(y)
            assert len(result) == len(y), "y and dydt are different lengths"
            
            for name in self.names_species:
                assert not np.any(np.isnan(self.y.current().get_species(name))), "NaN value in concentrations of %s" % name
                assert not np.any(np.isinf(self.y.current().get_species(name))), "Inf value in concentrations of %s" % name
            
            return result
            
        self.derivative = derivative
        
    def set_constant(self, name, value, split_axes=None):
        super(SimulationMultiprocessing, self).set_constant(name, value)
        self.variable_split_axes[name] = split_axes
        
    def set_variable_split_axes(self, dic):
        self.variable_split_axes = dic
        
    def split(self, name, value, default_slice):
        if name in sim.variable_split_axes.keys():
            slice_ = tuple([default_slice if x else slice(None) for x in sim.variable_split_axes[name]])
#            print value.shape, slice_
            res = value[slice_]
        else:
            res = value
        return res

def subset_slice(subset, superset):
    mask = np.in1d(superset, subset)
    ids = np.argwhere(mask)
    return slice(int(ids[0]), int(ids[-1]) + 1)
    
def init(dydt_rawarray_, ct_rawarray_, sim_):
    global dydt_rawarray, ct_rawarray, sim
    dydt_rawarray = dydt_rawarray_
    ct_rawarray = ct_rawarray_
    sim = sim_
        
def work(arg_vector):
    global dydt_rawarray, ct_rawarray, sim
    try:
        t, cids, inner_cids = arg_vector
        slice_cids = subset_slice(cids, sim.cids)
        slice_inner_cids = subset_slice(inner_cids, cids)
        
        ct_local = multicell.parallel.ConcentrationTableMultiprocessing(sim.names_species, cids, sim.cids, ct_rawarray)
        dydt_local = multicell.parallel.ConcentrationTableMultiprocessing(sim.names_species, inner_cids, sim.cids, dydt_rawarray)
        
        variables = dict((name, ct_local.get_species(name)) for name in sim.names_species)
        
        environment = variables
        environment["t"] = t
        environment["simulation"] = sim
        environment["np"] = np
        environment["math"] = math
        
        split_constants = {}
        for name, value in sim.parameters.items():
            split_constants[name] = sim.split(name, value, slice_cids)
        environment.update(split_constants)
        
        for name, expression in sim.intermediary_variables.items():
            environment[name] = eval(expression, environment)
            
        environment.update({"diffusion": multicell.simulation.diffusion,
                            "compute_differences": multicell.simulation.compute_differences,
                            "transport_against_gradient": multicell.simulation.transport_against_gradient,
                            "adjacency_matrix": sim.adjacency_matrix[slice_cids, slice_cids]})
        for name in sim.names_species:
            eval_result = eval(sim.derivative_components[name], environment)
            if not isinstance(eval_result, np.ndarray):
                eval_result = np.full(len(cids), eval_result, dtype=ctypes.c_double)
            elif eval_result.shape <> (len(cids),):
                eval_result *= np.ones(len(cids))
            dydt_local.set_species(name, eval_result[slice_inner_cids])
#    return eval_result
    except:
        raise Exception("".join(traceback.format_exception(*sys.exc_info())))
        
dydt_rawarray, ct_rawarray = None, None
sim = None