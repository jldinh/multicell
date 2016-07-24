# -*- coding: utf-8 -*-
"""
Created on Fri Oct 02 10:41:35 2015

@author: Jean-Louis Dinh
"""

import concentration_table

class Series_Concentration_Table(list):
    
    def __init__(self):
        super(Series_Concentration_Table, self).__init__()
        self.current_time_index = -1
    
    def append_new_table(self, variables_list, cids):
        self.append_table(concentration_table.Concentration_Table(variables_list, cids))
        
    def current(self):
        return self[self.current_time_index]
        
    def append_table(self, table):
        self.current_time_index += 1
        self.append(table)

    def get_species(self, name, cid=None, time_index=None):
        if time_index == None:
            time_index = self.current_time_index
        return self[time_index].get_species(name, cid)
    
    def set_species(self, name, values, time_index=None):
        if time_index == None:
            time_index = self.current_time_index
        self[time_index].set_species(name, values)
    

        