# -*- coding: utf-8 -*-
"""
Created on Sun Dec 26 17:59:27 2021

@author: Admin
"""

import numpy as np
import matplotlib.pyplot as plt


class Unit:
    def __init__ (self, type=None, length=1.0, Dax=0., ax_disc=1):
        self.type = type
        self.length = length
        self.Dax = Dax
        
        # Discretization
        self.ax_disc = ax_disc
        self.dz = self.length / self.ax_disc if self.ax_disc != 0. else 1.
        self.dz2 = self.dz ** 2
        
        self.c = np.zeros(self.ax_disc)
        

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
        self.time_section_index = 0
        self.system_parameters = None
        self.t = 0
        self.dt = 0.01
        self.c = np.zeros(50)
        
        
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
        ci = self.system_parameters.inlet_concentration
        
        self.c[0] = ci + (self.units[1].Dax / u) * (self.c[1] - self.c[0]) / self.units[1].dz
        
        self.c[1:-1] = self.c[1:-1] + \
            (self.units[1].Dax * ((self.c[2:] - 2*self.c[1:-1] + self.c[:-2]) / self.units[1].dz2)
                              -  u * ((self.c[1:-1] - self.c[:-2]) / self.units[1].dz)) * self.dt    
        
        self.c[-1] = self.c[-1] - (u * (self.c[-1] - self.c[-2]) / self.units[1].dz) * self.dt
        
        self.t += self.dt
    
    
    def run (self):
        total_time = sum([i[0] for i in self.time_sections])
        
        while self.t <= total_time:
            
            plt.clf()
            plt.title(f"t={round(self.t, 3)}")
            plt.ylim((0, 2))
            plt.plot(range(self.c.size), self.c)
            plt.pause(0.01)
            
            self.step()
            
        plt.show()