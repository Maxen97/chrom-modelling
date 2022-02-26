# -*- coding: utf-8 -*-
"""
Created on Sun Dec 26 17:59:27 2021

@author: Admin
"""

import numpy as np
import matplotlib.pyplot as plt

from unit import Unit


class Component:
    def __init__ (self):
        self.mw = 10000
        self.ads = 1
        self.kf = 1e-8
        self.v = 12
        self.s = 4
        
        
class Solution:
    def __init__ (self, components=[], concentrations=[]):
        self.components = components
        self.concentrations = concentrations
        

class Phase:
    def __init__ (self, inlet_solution, flowrate, length):
        self.inlet_solution = inlet_solution
        self.flowrate = flowrate
        self.length = length
    

class Experiment:
    def __init__(self):
        self.units = []
        self.t = 0
        self.dt = 0.01
        
        self.phases = []
        
        self.components = []
        
        # Phase-specific values
        self.f = 0
        self.components_inlet_concentrations = []
        
        self.t_solve = []
        self.c_solve = []
        
        
    def add_phase(self, phase):
        self.phases.append(phase)
    
    def update_phases(self):
        for phase in self.phases:
            if self.t <= phase.length:
                self.f = phase.flowrate
                self.components_inlet_concentrations = phase.inlet_solution.concentrations
                return
            
    
    def step (self):
        self.update_phases()
        
        c0 = self.components_inlet_concentrations
        
        c_previous = c0
        for unit in self.units:
            unit.step(c_previous, self.f, self.dt)
            c_previous = unit.c_out
            
        self.c_solve.append(list(self.units[0].c_out))
        
        self.t += self.dt
    
    
    def run (self):
        total_time = sum([i.length for i in self.phases])
        self.units[0].components = self.components
        
        while self.t <= total_time:
            self.step()
            self.t_solve.append(self.t)
        
        for i in range(len(self.components)):
            plt.plot(self.t_solve, [x[i] for x in self.c_solve])
            
        plt.show()