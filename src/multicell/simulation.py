# -*- coding: utf-8 -*-

import numpy as np
from odesparse import odeints
import scipy as sp
import re
from multicell import concentration_table, series_concentration_table
from openalea.tissueshape import centroid, divide_cell, cell_volume, face_surface_3D
import time
import sys
import math
#import numexpr as ne
from multicell.utilities import print_flush
#from multicell import utilities
import collections
import inspect
import functools

#ne.set_num_threads(8)

def restrict_environment(environment, function):
    parameters = inspect.getargspec(function)[0]
    return {k: environment[k] for k in parameters} 

# Compute a matrix of differences between cell concentrations for species name
def compute_differences(values):
    """
    Computes a matrix of pairwise differences between components of a vector.
    
    Parameters
    ----------
    values : ndarray
        1D vector containing the values to operate on.
        
    Returns
    -------
    ndarray
        Square ndarray containing the pairwise differences between components
        of value.
    """
    
    # horizontal -> column pattern = other
    # vertical vector -> row pattern = self
    # self - other
    return np.reshape(values, (-1, 1)) - values 

def diffusion(D, values, adjacency_matrix):
    """
    Computes diffusion in a tissue.
    
    The tissue is define by the given of its concentration values for the
    species of interest and a matrix of adjacency for its cells. Note that the
    matrix does not have to be the physical adjacency matrix, but can take into
    account the permeability (and impermeability) of cell interfaces (e.g.: if
    a given species can only diffuse in a subpopulation of cells).
    
    Parameters
    ----------
    D : float
        Diffusion coefficient
    values : ndarray
        1D vector containing the concentration values of the cells.
    adjacency_matrix : csr_matrix
        Square sparse matrix containing the exchange surfaces between any pair
        of cells in the tissue. Cells should be ordered in the same way as in
        values.
    
    Returns
    -------
    ndarray
        1D vector of chemical species flows, of the same length as values.
    """
    
    differences = compute_differences(values)
    return -D * (adjacency_matrix.multiply(differences)).sum(1).getA().flatten()

def transport_against_gradient(T, values, adjacency_matrix, response="linear", params={}):
    """
    Implementation of Jönsson, 2006 (simple version).
    
    The transporters of the species of interest are assumed to be in equal 
    quantities in all cells of the tissue. They relocate to cell interfaces
    proportionally to::
        * The surface of the interface;
        * The concentration in the species of interest on the other side of the
        interface.
    
    Parameters
    ----------
    T : float
        Transport coefficient.
    values : ndarray
        1D vector containing the concentrations of the cells of the tissue.
    adjacency_matrix : csr_matrix
        Square sparse matrix containing the exchange surfaces between any pair
        of cells in the tissue. Cells should be ordered in the same way as in
        values.
        
    Returns:
    --------
    ndarray
        1D vector of chemical species flows, of the same length as values.
    """
    
    # 1 value per column => horizontal broadcast into column pattern
    # Should be precomputed for incomplete grids
    sum_concentrations_neighbors = adjacency_matrix.multiply(values).sum(1).getA().flatten()
    sum_concentrations_neighbors += 1 * (sum_concentrations_neighbors == 0)
    
    inv_sum_concentrations_neighbors = 1. / sum_concentrations_neighbors
    
    incoming_transporters_distribution = adjacency_matrix.multiply(np.reshape(values, (-1, 1))).multiply(inv_sum_concentrations_neighbors)
    outgoing_transporters_distribution = adjacency_matrix.multiply(values).multiply(np.reshape(inv_sum_concentrations_neighbors, (-1, 1)))
    
    #assert abs(np.linalg.norm(transporters_distribution.sum(1), np.inf) - 1) < 0.001
    
    if response == "linear":
        values_response = values
    elif response == "hill":
        values_n = values ** params["n"]
        values_response = values_n / (params["threshold"] ** params["n"] + values_n)

    incoming_transport_fluxes = incoming_transporters_distribution * values_response
    outgoing_transport_fluxes = outgoing_transporters_distribution * np.reshape(values_response, (-1, 1))
    
    outgoing = outgoing_transport_fluxes.sum(1)
    incoming = incoming_transport_fluxes.sum(1)
    return T * (incoming - outgoing)

def parse_for_dependencies(functions, keys):
    """
    Finds which entities the expressions depends on.
    
    Parameters
    ----------
    expressions : dict(string -> string)
        Dictionary (entity -> expression) of expressions to parse for 
        dependencies.
    keys : list(string)
        List of entities whose functions should be analyzed.
        
    Returns
    -------
    dict(string, list(string))
        Dictionary (entity -> list of entities it depends on).
    """
    
    return dict(zip(keys, [inspect.getargspec(functions[k])[0] for k in keys]))

def null_function(*args, **kwargs):
    return 0.

class Simulation(object):
    """
    Class used to define all the information needed to run a tissue simulation.
    
    Contains generic methods to help with tissue simulations.
    """
    
    def __init__(self):
        """Initializes the Simulation with default values"""
        
        self.parameters = {}
        self.intermediary_variables = collections.OrderedDict()
        self.names_species = []
        self.intensive_cell_variables = set()
        self.n_species = 0
#        self.initialize_cell_grid()
        self.ODEs = {}
        self.variables_adjacency_matrices = {}
        self.adjacency_matrices = {}
        self.adjacency_matrices_specs = {"adjacency_matrix": lambda: self.cids}
        self.growth = False
        self.division = False
        self.render = False
        self.n_time_steps = 2
        self.n_steps_growth = 1
        self.volume_threshold_division = 1.5
        self.verbose = True
        self.dilution_volume_function = cell_volume
        self.dilution_volume_function_parameters = {}
        self.contraction = 0.
        self.t = [0., 1.]
        self.save_pictures = False
        self.detection_mode = False
        
    def set_verbose(self, verbose):
        """
        Turns verbose output on or off.
        
        Parameters
        ----------
        verbose : bool
            Whether to print the verbose output to the console or not.
        """
        
        self.verbose = verbose
    
    def import_zones(self, zones):
        """
        Imports zone information for the tissue.
        
        This information should follow the format defined in the tissuemeca
        module.
        
        Parameters
        ----------
        zones : dict(int, int)
            Dictionary mapping cell IDs to zone IDs.
        """
        
        self.zones = zones
    
    def set_parameter(self, name, value):
        """
        Deprecated. See set_constant.
        """
        
        self.set_constant(name, value)
    
    def set_constant(self, name, value):
        """
        Defines a constant with the given name and value.
        
        The constant can be used in ODE equations.
        
        Parameters
        ----------
        name : string
            Name of the constant.
        value : string
            Value of the constant.
        """
        
        self.parameters[name] = value
        
    def set_parameters(self, dictionary):
        """
        Deprecated. See set_constants.
        """
        
        self.set_constants(dictionary)
        
    def set_constants(self, dictionary):
        """
        Merge a dictionary of parameter values into that of the Simulation.
        
        Parameters
        ----------
        dictionary : dict(string, any)
            Dictionary containing pairs of names and values.
        """
        
        self.parameters.update(dictionary)
    
    def register_computed_variable(self, name, expression):
        """
        Defines a variable that is recomputed at each iteration.
        
        Computed variables defined this way are only computed once per
        iteration, but can be used multiple times in equations at no additional
        cost.
        
        Parameters
        ----------
        name : string
            Name of the computed variable.
        expression : string
            Python string used to compute the value of the variable.
        """
        
        self.intermediary_variables[name] = expression
        
    def make_dual_function(self, name, operator):
        code = "def dual_function(simulation, %s, **kwargs): return %s %s simulation.dilution_volumes" % (name, name, operator)
        exec(code)
        return locals()["dual_function"]
    
    def register_cell_variable(self, name, intensive=False, create_dual=True):
        """
        Register a variable of the given name that all cells will instantiate.
        
        Parameters
        ----------
        name : string
            Name of the variable to register
        intensive : bool
            If True, the variable is intensive, otherwise it is extensive. This
            affects how it is inherited during cell divisions.
        """
        
        if type(name) == list:
            for e in name:
                self.register_cell_variable(e, intensive)
        else:
            self.names_species.append(name)
            if intensive:
                self.intensive_cell_variables.add(name)
                if create_dual:
                    self.register_computed_variable("q_" + name, self.make_dual_function(name, "*"))
            else:
                if create_dual:
                    self.register_computed_variable("c_" + name, self.make_dual_function(name, "/"))
            self.set_ODE(name, null_function)
            self.n_species += 1
    
    def get_surface(self, wid):
        """
        Retrieves the surface of a face delimiting one or two cells.
        
        Parameters
        ----------
        wid : int
            ID of the interface in the Topomesh of the tissue.
        """
        
        return face_surface_3D(self.mesh, self.get_pos(), wid)
        
    def compute_dilution_volumes(self):
        """
        Computes the volumes used for concentration calculations from
        quantities of matter.
        
        Volumes are ordered like cell IDs.
        
        Returns
        -------
        PropertyTable
            Array of cell volumes.
        """
        
        self.dilution_volumes = concentration_table.Concentration_Table(["dilution_volumes"], self.cids)
        self.dilution_volumes.import_values(np.array([self.dilution_volume_function(self.mesh, self.get_pos(), cid, **self.dilution_volume_function_parameters) for cid in self.cids]))
        return self.dilution_volumes
        
    def compute_generic_adjacency_matrix(self, cids_function):
        """
        Compute the adjacency matrix (matrix of exchange surfaces) between the
        cells whose IDs are given.
        
        Parameters
        ----------
        cids : list or ndarray
            List of cells IDs for which the matrix should be computed.
            
        Returns
        -------
        lil_matrix
            Matrix containing, for all pairs of cells in the tissue, the
            exchange surface between the two cells if both are in cids, or 0
            otherwise.
        """
        
        cids = cids_function()
        adjacency_matrix = sp.sparse.lil_matrix((self.n_cells, self.n_cells), dtype=np.float64)
        for cid in cids:
            for wid in self.mesh.borders(3, cid):
                s = self.get_surface(wid)
                for ocid in self.mesh.regions(2, wid):
                    if ocid <> cid and ocid in cids:
                        adjacency_matrix[self.dict_cids[cid], self.dict_cids[ocid]] += s
        return adjacency_matrix
    
    def compute_adjacency_matrix(self):
        """
        Compute the adjacency matrix (matrix of exchange surfaces) between all
        cells of the tissue.
        
        Returns
        -------
        lil_matrix
            Sparse matrix of exchange surfaces
        """
        
        print_flush("Adjacency matrix computation: started")
        time_start = time.time()
        self.adjacency_matrices["adjacency_matrix"] = self.compute_generic_adjacency_matrix(self.cids).tocsr()
        print_flush("Adjacency matrix computation: finished (%.2f s)" % (time.time() - time_start))
    
    def register_adjacency_matrix(self, name, cids_function):
        self.adjacency_matrices_specs[name] = cids_function
    
    def compute_adjacency_matrices(self):
        for name, cids_function in self.adjacency_matrices_specs.items():
            self.adjacency_matrices[name] = self.compute_generic_adjacency_matrix(cids_function)
    
    def _transport(self, function, coeff, values, adjacency_matrix, *args, **kwargs):
        if self.detection_mode:
            assert len(values.variables_list) == 1
            name = values.variables_list[0]
            self.add_adjacency_matrix(name, adjacency_matrix)
        return function(coeff, values, adjacency_matrix, *args, **kwargs)
    
    def diffusion(self, D, values, adjacency_matrix):
        return self._transport(diffusion, D, values, adjacency_matrix)
    
    def transport_against_gradient(self, T, values, adjacency_matrix, *args, **kwargs):
        return self._transport(transport_against_gradient, T, values, adjacency_matrix, *args, **kwargs)
    
    def add_adjacency_matrix(self, name, adjacency_matrix):
        if name in self.variables_adjacency_matrices:
            self.variables_adjacency_matrices[name] += adjacency_matrix
        else:
            self.variables_adjacency_matrices[name] = adjacency_matrix
    
    def compute_Jacobian(self):
        """
        Compute the Jacobian matrix of the ODE system.
        """
        
        print_flush("Jacobian computation: started")
        if not hasattr(self, "dependencies"):
            self.compute_dependencies()
        time_start = time.time()
        #self.compute_adjacency_matrix()
        n_cells = self.n_cells
        n_species = self.n_species
        size_JPat = n_species*n_cells
        JPat = sp.sparse.lil_matrix((size_JPat, size_JPat), dtype=np.float64)
        
        # Transport
        for i, v in enumerate(self.names_species):
            if v in self.variables_adjacency_matrices:
                matrix = self.variables_adjacency_matrices[v]
                JPat[i*n_cells:(i+1)*n_cells, i*n_cells:(i+1)*n_cells] = matrix
            
        # Cell-autonomous effects
        for i in xrange(n_species):
            for k in xrange(n_cells):
                for j in [self.names_species.index(effector) for effector in self.dependencies[self.names_species[i]] if effector in self.names_species]:
                    JPat[i*n_cells+k, j*n_cells+k] = 1.
                    
        self.JPat = JPat
        
        print_flush("Jacobian computation: finished (%.2f s)" % (time.time() - time_start))
    
    def set_mesh(self, mesh):
        """
        Sets the given mesh as the mesh used in the simulation.
        
        
        Parameters
        ----------
        mesh : Topomesh
            Topomesh of the tissue to use in the simulation.
        """
        
        self.mesh = mesh
        
    def initialize_mesh_properties(self):
        """
        Compute convenient properties related to the Topomesh of the tissue.
        
        These properties include the IDs of cells, the number of cells and a
        dictionary mapping cell IDs to their position in the various vectors
        and matrices.
        They are stored in the Simulation object as members (cids, n_cells and
        dict_cids, respectively)
        """
        
        self.cids = np.array(list(self.mesh.wisps(3)))
        self.n_cells = len(self.cids)
        self.dict_cids = dict((cid, i) for (i, cid) in enumerate(self.cids))
        self.compute_dilution_volumes()
    
    def get_pos(self):
        """
        Returns the positions of the vertices of the tissue.
        
        Returns
        -------
        dict(int, ndarray)
            Dictionary mapping vertex IDs to arrays of coordinates (1 array per
            ID).
        """
        
        return self.pos
    
    def set_pos(self, pos):
        """
        Sets the positions of the vertices of the tissue.
        
        Parameters
        ----------
        pos : dict(int, ndarray)
            Dictionary mapping vertex IDs to arrays of coordinates (1 array per
            ID).
        """
        
        self.pos = pos
    
    def import_topomesh(self, mesh, pos):
        """
        Import a Topomesh and the positions of its vertices
        
        Parameters
        ----------
            mesh : Topomesh
                Topomesh of the tissue to import into the Simulation.
            pos : dict(int, ndarray)
                Dictionary mapping vertex IDs to their positions.
        """
        
        time_start = time.time()
        print_flush("Topomesh importation: started")
        print_flush("- setting mesh")
        self.set_mesh(mesh)
        print_flush("- setting pos")
        self.set_pos(pos)
        print_flush("- updating properties")
        self.initialize_mesh_properties()
        print_flush("Topomesh importation: finished (%.2f s)" % (time.time() - time_start))
    
    def set_duration(self, t_max):
        """
        Sets the duration that the simulation should run for.
        
        Parameters
        ----------
            t_max : float
                Duration of the simulation.
        """
        
        self.t = [0., float(t_max)]
    
    def set_cell_variable(self, name, values, time_index=None):
        self.y.set_species(name, values, time_index)
        
    def get_cell_variable(self, name, cid=None, time_index=None):
        return self.y.get_species(name, cid, time_index)
    
    def initialize_concentrations(self, name=None, _min=0., _max=0.):
        """
        Deprecated.
        """
        
        if not hasattr(self, "y") or len(self.y) < 1 or self.y.current().shape <> (self.n_species, self.n_cells):
            self.y = series_concentration_table.Series_Concentration_Table()
            self.y.append_new_table(self.names_species, self.cids)
            
        if name == None:
            self.y[0][...] = (np.random.uniform(_min, _max, self.n_cells * self.n_species)
            .reshape((self.n_species, self.n_cells)))
        else:
            self.y.set_species(name, np.random.uniform(_min, _max, self.n_cells), 0)
    
    def initialize_cell_variables(self):
        """
        Creates the time series of cell variables and initializes all values to
        0.
        """
        
        if not hasattr(self, "y") or len(self.y) < 1 or self.y.current().shape <> (self.n_species, self.n_cells):
            self.y = series_concentration_table.Series_Concentration_Table()
            self.y.append_new_table(self.names_species, self.cids)
    
    def set_ODE(self, name, function):
        """
        Defines the ODE used to compute a cell variable.
        
        Parameters
        ----------
        name : string
            Name of the variable.
        function : function
            ODE of the variable describing cell-autonomous regulations.
        adjacency_matrices : [str, ...]
            List of the names assigned to the relevant adjacency matrices.
        """
        
        self.ODEs[name] = function
    
    def compute_dependencies(self):
        """
        Computes the dependencies of each ODE in the system.
        
        Returns
        -------
        dict(str, set)
            Dictionary mapping the names of cell variables to the other
            variables it depends on.
        """
        
        dependencies = parse_for_dependencies(self.ODEs, self.names_species)
        indirect_dependencies = parse_for_dependencies(self.intermediary_variables, self.intermediary_variables.keys())
        new_dependencies = {}
        
        while new_dependencies <> dependencies:
            new_dependencies = dependencies
            for target in self.names_species:
                for (i_regulator, regulator) in enumerate(dependencies[target]):
                    if regulator in self.intermediary_variables.keys():
                        new_dependencies[target].pop(i_regulator)
                        new_dependencies[target] += indirect_dependencies[regulator]
        
        #Filtering
        for name in self.names_species:
            dependencies[name] = set(dependencies[name]) & set(self.names_species)
        
        self.dependencies = dependencies
    
    def compute_environment(self, y, t):
        ct = concentration_table.Concentration_Table(self.names_species, self.cids).import_values(y)
        ct = ct * (ct>0)
        variables = dict((name, ct.get_species(name)) for name in self.names_species)
            
        environment = variables
        environment["t"] = t
        environment["simulation"] = self
        environment.update(self.parameters)
        for name, function in self.intermediary_variables.items():
            environment[name] = function(**restrict_environment(environment, function))
        environment.update(self.adjacency_matrices)
        self.derivative_environment = environment
        
    def detect_adjacency_matrices(self):
        self.detection_mode = True
        self.variables_adjacency_matrices = {}
        self.test_derivative()
        self.detection_mode = False
        
    def prepare_derivative(self):
        """
        Compiles the expressions of ODEs and prepares a function to feed into
        the solver.
        """
        
#        self.compute_dilution_volumes()
        self.compute_adjacency_matrices()
            
        def derivative(y, t):
            # Initialization of the derivative vector
            dydt = concentration_table.Concentration_Table(self.names_species, self.cids)
            
            self.compute_environment(y, t)
            for name in self.names_species:
                function = self.ODEs[name]
                dydt.set_species(name, function(**restrict_environment(self.derivative_environment, function)))
                
            result = dydt.as_1d_array()

            # Test
            assert len(result) == len(y), "y and dydt are different lengths"
            
            for name in self.names_species:
                assert not np.any(np.isnan(self.y.current().get_species(name))), "NaN value in concentrations of %s" % name
                assert not np.any(np.isinf(self.y.current().get_species(name))), "Inf value in concentrations of %s" % name
            
            return result
            
        self.derivative = derivative
        
        self.compute_dependencies()
        self.detect_adjacency_matrices()
        self.compute_Jacobian()
        
    def simulate(self):
        """
        Runs the simulation.
        """
        
#        self.prepare_derivative()
#        n_cells = self.n_cells
#        ts_initial_Jacobian = time.time()
#        self.compute_Jacobian()
#        print_flush("Initial computation of the Jacobian: %s seconds" % (time.time() - ts_initial_Jacobian))
        
        t = np.linspace(self.t[0], self.t[1], num=self.n_steps_growth+1)
        for i in xrange(self.n_steps_growth):
            self.prepare_derivative()
            if self.growth:
                print_flush("Growth step #%s" % i)
                
            ts_integration = time.time()
            lrw = 100000000
            (output, self.odeints_debug) = odeints(self.derivative, self.y.current().as_1d_array(), t[i:i+2], mxstep=10000, lrw=lrw, JPat=self.JPat, full_output=1)
            #output = odeints(self.derivative, self.y.current().as_1d_array(), t[i:i+2], nnz=self.JPat.getnnz(), lrw=1000000, JPat=self.JPat)
            
            lrw_needed = self.odeints_debug["leniw"] + self.odeints_debug["lenrw"]
            while lrw_needed > lrw:
                print_flush("Work array of insufficient length for successful integration. Increasing array length to %s." % lrw)
                lrw = lrw_needed
                (output, self.odeints_debug) = odeints(self.derivative, self.y.current().as_1d_array(), t[i:i+2], lrw=lrw, JPat=self.JPat, full_output=1)
                lrw_needed = self.odeints_debug["leniw"] + self.odeints_debug["lenrw"]
                
            result = output[-1]
            self.y.append_table(concentration_table.Concentration_Table(self.names_species, self.cids).import_values(result))
            print_flush("Integration of the ODE system: %s seconds" % (time.time() - ts_integration))

            if self.growth:
                ts_growth = time.time()
                self.set_pos(self.growth_method(self.mesh, self.get_pos(), **self.growth_method_parameters))
                print_flush("Growth of the tissue: %s seconds" % (time.time() - ts_growth))
                
            if self.division:
                ts_division = time.time()
                list_divisions = self.divide_all()
                #self.update_Jacobian(list_divisions)
                print_flush("Cell divisions: %s seconds" % (time.time() - ts_division))
            
#            if self.growth or self.division:
#                self.compute_Jacobian()
#                self.compute_dilution_volumes()
            
            if self.render:
                self.renderer.display(self.rendered_species, save=self.save_pictures)
    
    def enable_growth(self, n_steps):
        """
        Enables growth.
        """
        
        self.growth = True
        self.n_steps_growth = n_steps
    
    def enable_division(self, contraction=0.):
        """
        Enables cell divisions.
        """
        
        self.division = True
        self.contraction = contraction
    
    def disable_growth(self):
        """
        Disables growth.
        """
        
        self.growth = False
        
    def disable_division(self):
        """
        Disables cell divisions.
        """
        
        self.division = False
    
    def register_growth_method(self, growth_method, growth_method_parameters={}):
        """
        Registers a growth method to be used during the simulation.
        
        Parameters
        ----------
        growth_method : function
            Function that makes a mesh, a positions dictionary and, optionally,
            other arguments, and returns a new pos dictionary.
        growth_method_parameters: dict(string, any)
            Dictionary of additional parameters required by growth_method
        """
        
        self.growth_method = growth_method
        self.growth_method_parameters = growth_method_parameters
    
    def register_division_method(self, division_method, division_method_parameters={}):
        """
        Registers a division method to be used during the simulation.

        Parameters
        ----------
        division_method : function
            Function that takes a mesh, a positions dictionary, a cell ID,
            and, optionally,
            other arguments, and returns a point in the division plane and a
            vector normal to the division plane.
        division_method_parameters : dict(string, any)
            Dictionary of additional parameters required by division_method
        """
        
        self.division_method = division_method
        self.division_method_parameters = division_method_parameters
        
    def register_division_trigger(self, division_trigger, division_trigger_parameters={}):
        """
        Registers a division trigger to be used during the simulation.

        Parameters
        ----------
        division_trigger : function
            Function that takes a mesh, a positions dictionary, a cell ID,
            and, optionally,
            other arguments, and returns a point in the division plane and a
            vector normal to the division plane.
        division_trigger_parameters : dict(string, any)
            Dictionary of additional parameters required by division_method
        """
        
        self.division_trigger = division_trigger
        self.division_trigger_parameters = division_trigger_parameters
    
    def register_dilution_volume_function(self, dilution_volume_function, dilution_volume_function_parameters={}):
        """
        Registers a volume function to be used during the simulation to compute
        species concentrations, and recomputes the volumes.

        Parameters
        ----------
        dilution_volume_function : function
            Function that takes a mesh, a positions dictionary, a cell ID,
            and, optionally,
            other arguments, and returns the volume to use to compute the
            concentrations of species in that cell.
        dilution_volume_function_parameters : dict(string, any)
            Dictionary of additional parameters required by division_method
        """
        
        self.dilution_volume_function = dilution_volume_function
        self.dilution_volume_function_parameters = dilution_volume_function_parameters
#        self.dilution_volumes = self.compute_dilution_volumes()
    
    #Adapted from http://openalea.gforge.inria.fr/doc/vplants/tissue/doc/_build/html/user/cell_division/simple3D/index.html
    def is_internal_point(self, pid):
        """
        Finds out whether a vertex is internal or not.
        
        Parameters
        ----------
        pid : int
            ID of the vertex.
        
        Returns
        -------
        bool
            Whether the vertex is internal or not.
        """
        
        for eid in self.mesh.regions(0, pid) :
            for fid in self.mesh.regions(1,eid) :
                if self.mesh.nb_regions(2,fid) == 1 :
                    return False
    	return True
    
    def is_external_edge(self, eid):
        """
        Finds out whether a edge is internal or not.
        
        Parameters
        ----------
        eid : int
            ID of the edge.
        
        Returns
        -------
        bool
            Whether the edge is internal or not.
        """
        
        for fid in self.mesh.regions(1,eid) :
            if self.mesh.nb_regions(2,fid) == 1 :
                return True
        return False
    
    def cell_division(self, cid, point, normal):
        """
        Performs the division of the specified cell, by a plane defined by a
        point and a normal. Optionally, shrinks the new cell wall with respect
        to the cell section, so the daughter cells look more natural.
        
        Parameters
        ----------
        cid : int
            ID of the cell to divide.
        point : ndarray
            Coordinates of a point in the division plane.
        normal : ndarray
            Coordinates of a normal vector of the division plane.
        shrink : float
            Coefficient between 0 and 1 indicating how much the new wall should
            shrink.
            
        Returns
        -------
        int
            ID of the mother cell.
        tuple(int, int)
            IDs of the daughter cells.
        """
        
        #perform division
        lineage = divide_cell(self.mesh, self.get_pos(), cid, point, normal, 0.)
        
        if self.contraction > 0.:
            shrink = self.contraction
            #shrink newly formed separation wall
            sca = (100 - shrink) / 100.
            wid, = lineage[2][None]
            bary = centroid(self.mesh, self.get_pos(), 2, wid)
            for pid in self.mesh.borders(2,wid,2) :
                if self.is_internal_point(pid) :
                    self.get_pos()[pid] = bary + (self.get_pos()[pid] - bary) * sca
                    
            for eid in self.mesh.borders(2,wid) :
                if self.is_external_edge(eid) :
                    bary = centroid(self.mesh, self.get_pos(), 1, eid)
                    for pid in self.mesh.borders(1,eid) :
                        self.get_pos()[pid] = bary + (self.get_pos()[pid] - bary) * sca
                    
        #update properties
        return (cid, lineage[3][cid])
        
    def divide_all(self):
        """
        Checks the division trigger on all cells and divide those that meet the
        requirements.
        
        Returns
        -------
        list(tuple(int, tuple(int, int)))
            List of divisions that have occurred.
        """
        
        list_divisions = []
        for cid in np.random.permutation(self.cids):
            #if cell_volume(self.mesh, self.get_pos(), cid) > self.volume_threshold_division:
            if self.division_trigger(self.mesh, self.get_pos(), cid, **self.division_trigger_parameters):
                bary, axis = self.division_method(self.mesh, self.get_pos(), cid, **self.division_method_parameters)
                list_divisions.append(self.cell_division(cid, bary, axis))
        if list_divisions <> []:
            position = dict((cid, i) for (i, cid) in enumerate(self.cids))
            list_divisions = sorted(list_divisions, key=lambda x: position[x[0]])
            self.update_properties_after_divisions(list_divisions)
        return list_divisions
    
    def update_properties_after_divisions(self, list_divisions):
        added_cids = np.array(sum([new_cids for (cid, new_cids) in list_divisions], ()))
        removed_cids = np.array([cid for (cid, new_cids) in list_divisions])
        
        assert len(added_cids) == 2 * len(removed_cids)
        
        current_cids = self.cids
        conserved_cids = np.array([cid for cid in current_cids if cid not in removed_cids])
        updated_cids = np.hstack((current_cids[np.in1d(current_cids, removed_cids, invert=True)], added_cids))
        self.set_cids(updated_cids)
        variables_list = self.names_species
        
        volumes = self.compute_dilution_volumes()
        new_cells_volumes = volumes.restrict(added_cids).as_1d_array()
        grouped_volumes = np.repeat(np.sum(new_cells_volumes.reshape((2, -1), order='F'), axis=0), 2)
        
        extensive_coeffs = new_cells_volumes / grouped_volumes
        intensive_coeffs = np.ones(extensive_coeffs.shape)
        coeffs = np.vstack([intensive_coeffs if x in self.intensive_cell_variables else extensive_coeffs for x in self.names_species])
        new_cells_properties = np.repeat(self.y.current().restrict(removed_cids), 2, axis=1) * coeffs
        
        updated_values = np.hstack((self.y.current().restrict(conserved_cids), new_cells_properties))
        self.y.append_new_table(variables_list, updated_cids)
        self.y.current().import_values(updated_values)
    
    def set_cids(self, cids):
        self.cids = np.array(cids)
        self.n_cells = len(cids)
        self.dict_cids = dict((cid, i) for (i, cid) in enumerate(cids))
    
    def register_renderer(self, renderer_class, rendered_species=None):
        self.renderer = renderer_class(self)
        self.render = True
        self.rendered_species = rendered_species
    
    def test_derivative(self):
        self.derivative(self.y.current().as_1d_array(), 0)
