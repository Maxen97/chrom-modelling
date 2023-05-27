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
        self.cl = np.ones((len(self.components), self.ax_disc), dtype=np.longfloat)
        self.cp = np.ones((len(self.components), self.ax_disc), dtype=np.longfloat)
        self.q = np.ones((len(self.components), self.ax_disc), dtype=np.longfloat)
        
        # Concentrations at unit inlet
        self.c_in = np.zeros(len(self.components))
        
        # Concentrations at unit outlet
        self.c_out = np.zeros(len(self.components))
        
    def set_init_concentrations (self, c):
        conc = np.array(c)
        
        self.cl *= conc[:, None]
        self.cp *= conc[:, None]
        self.q *=  conc[:, None]
        
        
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
         
        """
        # Update bead concentrations
        for i in range(len(self.components)):
            self.cp[i][:] = self.cp[i][:] + \
                (3 * self.components[i].kf / (self.ep * self.rp * 1) * (self.cl[i][:] - self.cp[i][:])) * dt
        """
        
        # Update adsorption
        
        # Calculate number of bound counter-ions
        q0 = self.ic - sum([self.components[i].v * self.q[i] for i in range(len(self.components[1:]))])
        
        # Calculate number of free binding sites
        q0_hat = q0 - sum([self.components[i].s * self.q[i] for i in range(len(self.components[1:]))])
        
        qref = 1
        cref = 1
        
        for i in range(len(self.components) - 1):
            i =+ 1
            
            self.q[i, :] = self.q[i, :] + ((self.components[i].ka * self.cp[i, :] * (q0/qref)**self.components[i].v) - (self.components[i].kd * self.q[i, :] * (self.cp[0]/cref)**self.components[i].v)) * dt
            #sum_sigma = self.component
        
        
        # Update bead concentrations
        for i in range(len(self.components)):
            self.cp[i][:] = self.cp[i][:] + \
                (3 * self.components[i].kf / (self.ep * self.rp * 1) * (self.cl[i][:] - self.cp[i][:]) * \
                 (self.components[i].ka * self.cp[i, :] * (q0/qref)**self.components[i].v) - (self.components[i].kd * self.q[i, :] * (self.cp[0]/cref)**self.components[i].v)) * dt