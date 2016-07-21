# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 12:06:18 2015

@author: jl
"""

from multicell.simulation import Simulation
import time

class Simulation_MPI(Simulation):
    
    def __init__(self):
        super(Simulation_MPI, self).__init__()
        
    def set_intra_comm(self, intra_comm):
        self.intra_comm = intra_comm
        self.mpi_size = intra_comm.Get_size()
        self.mpi_rank = intra_comm.Get_rank()
        
    def set_inter_comm(self, inter_comm):
        self.inter_comm = inter_comm
        
    def listen_on_inter_comm(self):
        inter_comm = self.inter_comm
        while not inter_comm.Iprobe(0, 0):
            time.sleep(0.01)    
        inter_comm.recv(None, 0, 0)
        print str(self.mpi_rank) + " received start signal from Sofameca"