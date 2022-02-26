# -*- coding: utf-8 -*-
"""
Created on Fri Jan  7 15:12:41 2022

@author: Max
"""

import numpy as np


class Unit:
    def __init__ (self, type="NONE", length=0., Dax=0, ax_disc=1, components=[]):
        self.type = type
        self.length = length
        self.r = 0.010
        self.Dax = Dax
        self.rp = 20 * 1e-6
        self.et = 0.86
        self.ec = 0.37
        self.ep = (self.et - self.ec) / (1 - self.ec)
        
        self.bc = self.ec / (1 - self.ec)
        self.bp = self.ep / (1 - self.ep)
        
        
        self.ic = 140
        
        self.components = components
        
        # Discretization
        self.ax_disc = ax_disc
        self.dz = self.length / (self.ax_disc + 1)
        self.dz2 = self.dz ** 2
        
        # Concentration of each component at each node
        self.cl = np.zeros((len(self.components), self.ax_disc))
        self.cp = np.zeros((len(self.components), self.ax_disc))
        self.q = np.zeros((len(self.components), self.ax_disc))
        
        # Concentrations at unit inlet
        self.c_in = np.zeros(len(self.components))
        
        # Concentrations at unit outlet
        self.c_out = np.zeros(len(self.components))
        
    def set_init_concentrations (self, cl, cp, q):
        
        self.cl = np.zeros((len(self.components), self.ax_disc))
        self.cp = np.zeros((len(self.components), self.ax_disc))
        self.q = np.zeros((len(self.components), self.ax_disc))
        
        
    def step (self, c0, f, dt):
        # Convert from volumetric flow to linear interstitial flow
        u = f / (3.1415 * self.r**2 * self.ec)
        
        # Update bulk concentrations
        for i in range(len(self.components)):
            self.c_in = c0
            
            self.cl[i, 0] = self.c_in[i] + ((self.Dax / u) * (self.cl[i, 1] - self.cl[i, 0]) / (self.dz * 0.5)) * dt
            
            self.cl[i, 1:-1] = self.cl[i, 1:-1] + \
                (self.Dax * ((self.cl[i, 2:] - 2*self.cl[i, 1:-1] + self.cl[i, :-2]) / self.dz2)
                                  -  u * ((self.cl[i, 1:-1] - self.cl[i, :-2]) / self.dz) - \
                    3 * self.components[i].kf / (self.bc * self.rp) * (self.cl[i, 1:-1] - self.cp[i, 1:-1])) * dt
                    
            self.cl[i, -1] = self.cl[i, -1] - (u * (self.cl[i, -1] - self.cl[i, -2]) / (self.dz * 0.5)) * dt
                   
            self.c_out[i] = self.cl[i, -1]
         
        
        # Update bead concentrations
        for i in range(len(self.components)):
            self.cp[i][:] = self.cp[i][:] + \
                (3 * self.components[i].kf / (self.ep * self.rp * 1) * (self.cl[i][:] - self.cp[i][:])) * dt
        
        
        # Update adsorption
        sum_sigma = 0 # move inside loop
        for i in range(len(self.components)):
            pass
            #sum_sigma = self.component
        
        #q0 = self.ic - sum()