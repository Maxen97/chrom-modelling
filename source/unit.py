# -*- coding: utf-8 -*-
"""
Created on Fri Jan  7 15:12:41 2022

@author: Max
"""

import numpy as np


class Unit:
    def __init__ (self, type="NONE", length=0., Dax=0, ax_disc=1):
        self.type = type
        self.length = length
        self.Dax = Dax
        
        # Discretization
        self.ax_disc = ax_disc
        self.dz = self.length / self.ax_disc if self.ax_disc != 0. else 1.
        self.dz2 = self.dz ** 2
        
        # Concentration at nodes
        self.c = np.zeros(self.ax_disc)
        
        # Concentration at unit inlet
        self.c_in = 0
        
        # Concentration at unit outlet
        self.c_out = 0
        
        #test
    def step (self, c0, u, dt):
        self.c_in = c0
        
        # Left boundary condition
        self.c[0] = self.c_in + (self.Dax / u) * (self.c[1] - self.c[0]) / self.dz
        
        self.c[1:-1] = self.c[1:-1] + \
            (self.Dax * ((self.c[2:] - 2*self.c[1:-1] + self.c[:-2]) / self.dz2)
                              -  u * ((self.c[1:-1] - self.c[:-2]) / self.dz)) * dt    
        
        # Right boundary condition
        self.c[-1] = self.c[-1] - (u * (self.c[-1] - self.c[-2]) / self.dz) * dt
        
        self.c_out = self.c[-1]