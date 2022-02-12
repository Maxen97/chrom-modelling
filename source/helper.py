# -*- coding: utf-8 -*-
"""
Created on Sun Dec 26 17:59:27 2021

@author: Admin
"""

import numpy as np
import matplotlib.pyplot as plt

from unit import Unit


class TimeSection:
    def __init__ (self, system_parameters=None, time=0):
        self.time = time
        self.system_parameters = system_parameters
    

class SystemParameters():
    def __init__ (self):
        self.u = 0
        self.inlet_concentration = 0
    

class Experiment:
    def __init__(self):
        self.model = None
        self.units = []
        self.time_sections = []
        self.system_parameters = None
        self.t = 0
        self.dt = 0.001
        self.c = np.zeros(50)
        self.cp = np.zeros(50)
        
        self.components = []
        self.components_inlet_concentrations = []
        
        self.t_solve = []
        self.c_solve = []
        
    def add_timesection (self, time_section):
        self.time_sections.append((time_section.time, time_section.system_parameters))
        
    
    def update_timesection(self):
        for i in self.time_sections:
            if self.t <= i[0]:
                self.system_parameters = i[1]
                return
            
    
    def step (self):
        self.update_timesection()
        
        u = self.system_parameters.u
        c0 = self.system_parameters.inlet_concentration
        
        c_previous = c0
        for unit in self.units:
            unit.step(c_previous, u, self.dt)
            c_previous = unit.c_out
        self.c = self.units[0].cl
        
        self.t += self.dt
    
    
    def run (self):
        total_time = sum([i[0] for i in self.time_sections])
        
        while self.t <= total_time:
            self.step()
            self.t_solve.append(self.t)
            self.c_solve.append(self.c[-1])
        
        plt.plot(self.t_solve, self.c_solve)
        plt.show()