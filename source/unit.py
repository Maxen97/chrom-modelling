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
        self.r = 0.010
        self.Dax = Dax
        self.rp = 20 * 1e-6
        self.et = 0.85
        self.ec = 0.37
        self.ep = (self.et - self.ec) / (1 - self.ec)
        
        self.bc = self.ec / (1 - self.ec)
        
        self.kf = 1e-8
        
        # Discretization
        self.ax_disc = ax_disc
        self.dz = self.length / self.ax_disc if self.ax_disc != 0. else 1.
        self.dz2 = self.dz ** 2
        
        # Concentration at nodes
        self.cl = np.zeros(self.ax_disc)
        self.cp = np.zeros(self.ax_disc)
        
        # Concentration at unit inlet
        self.c_in = 0
        
        # Concentration at unit outlet
        self.c_out = 0
        
        
    def step (self, c0, f, dt):
        u = f / (3.1415 * self.r**2 * self.ec)
        
        self.c_in = c0
        
        # Left boundary condition
        self.cl[0] = self.c_in + (self.Dax / u) * (self.cl[1] - self.cl[0]) / self.dz
        
        self.cl[1:-1] = self.cl[1:-1] + \
            (self.Dax * ((self.cl[2:] - 2*self.cl[1:-1] + self.cl[:-2]) / self.dz2)
                              -  u * ((self.cl[1:-1] - self.cl[:-2]) / self.dz)) * dt - \
                1 / self.bc * 3 * self.kf / self.rp * (self.cl[1:-1] - self.cp[1:-1])
        
        # Right boundary condition
        self.cl[-1] = self.cl[-1] - (u * (self.cl[-1] - self.cl[-2]) / self.dz) * dt
        
        self.c_out = self.cl[-1]
        
        
        # Update cp
        self.cp[1:-1] = self.cp[1:-1] + 3 * self.kf / (self.ep * self.rp) * (self.cl[1:-1] - self.cp[1:-1])