# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 20:50:50 2015

@author: jl
"""

import numpy as np

class Concentration_Table(np.ndarray):
    
    def __new__(cls, variables_list, cids):
        zeros = np.zeros((len(variables_list), len(cids)), dtype=np.float64)
        array = np.ndarray.__new__(cls, zeros.shape, zeros.dtype, np.getbuffer(zeros))
        array.variables_list = variables_list
        array.cids_array = np.array(cids)
        array._rebuild_cached_info()
        return array
        
    def __array_finalize__(self, a):
        if a is None:
            return
        self.variables_list = getattr(a, 'variables_list', None)
        self.cids_array = getattr(a, 'cids_array', None)
        self._rebuild_cached_info()
        
    def __getnewargs__(self):
        return (self.variables_list, self.cids_array)
        
    def __copy__(self):
        return self
        
    def __deepcopy__(self, memo):
        return self
        
    def __reduce__(self):
        # Get the parent's __reduce__ tuple
        pickled_state = super(Concentration_Table, self).__reduce__()
        # Create our own tuple to pass to __setstate__
        new_state = pickled_state[2] + (self.variables_list, self.variables_dict, self.cids_array, self.dict_cids, self.n_cells)
        # Return a tuple that replaces the parent's __setstate__ tuple with our own
        return (pickled_state[0], pickled_state[1], new_state)

    def __setstate__(self, state):
        specific_state = state[-5:]
        self.variables_list, self.variables_dict, self.cids_array, self.dict_cids, self.n_cells = specific_state  # Set the info attribute
        # Call the parent's __setstate__ with the other tuple elements.
        super(Concentration_Table, self).__setstate__(state[0:-5])
        
    def _rebuild_cached_info(self):
        self.variables_dict = dict((name, i) for (i, name) in enumerate(self.variables_list))
        self.dict_cids = dict((cid, i) for (i, cid) in enumerate(self.cids_array))
        self.n_cells = len(self.cids_array)
        
    def get_species(self, name, cid=None):
        if cid == None:
            res = self[self.variables_dict[name], :]
            res.variables_list = [name]
            res.variables_dict = {name: 0}
            res._rebuild_cached_info()
            return res
        else:
            return self.get_species(name)[self.dict_cids[cid]]
    
    def set_species(self, name, values):
        self[self.variables_dict[name], :] = values
    
    def as_1d_array(self):
        dim = np.prod(self.shape)
        return np.reshape(self, (dim,), 'C')

    def import_values(self, values):
        self[...] = values.reshape(self.shape)
        return self
    
    def get_indices(self, cids):
        return np.array([self.dict_cids[cid] for cid in cids], dtype=np.int32)
    
    def restrict(self, cids):
        indices = self.get_indices(cids)
        return self[:, indices]