import numpy as np
import matplotlib.pyplot as plt

class LRMP:
    """ Lumped rate model with pores
    """
    def __init__ (self, nz, lz, lr, et, ec, rp, dax, components, initial_concentrations):
        self.nz = nz
        self.lz = lz
        self.lr = lr
        self.et = et
        self.ec = ec
        self.rp = rp
        self.dax = dax
        self.components = components
        self.initial_concentrations = initial_concentrations

        self.ep = (self.et - self.ec) / (1 - self.ec)
        self.bc = self.ec / (1 - self.ec)
        self.bp = self.ep / (1 - self.ep)

        # Discretization
        self.dz = self.lz / (self.nz + 1)
        self.dz2 = self.dz ** 2
        
        # Concentration matrices
        self.cl = np.zeros((len(self.components), self.nz), dtype=np.longfloat) # interstitial liquid
        self.cp = np.zeros((len(self.components), self.nz), dtype=np.longfloat) # pore liquid
        self.c_in = np.zeros(len(self.components), dtype=np.longfloat) # unit inlet
        self.c_out = np.zeros(len(self.components), dtype=np.longfloat) # unit outlet
        
    def set_init_concentrations (self, c):
        conc = np.array(c)
        
        self.cl *= conc[:, None]
        self.cp *= conc[:, None]
        self.q *=  conc[:, None]


    def step (self, c0, f, dt):
        # Convert from volumetric flow to interstitial linear flow
        u = f / (np.pi * self.lr**2 * self.ec)
        
        # Update bulk concentrations
        for i in range(len(self.components)):
            self.c_in = c0
            
            self.cl[i, 0] = self.c_in[i] + ((self.dax / u) * (self.cl[i, 1] - self.cl[i, 0]) / (self.dz * 0.5)) * dt
            
            self.cl[i, 1:-1] = self.cl[i, 1:-1] + \
                (self.dax * ((self.cl[i, 2:] - 2*self.cl[i, 1:-1] + self.cl[i, :-2]) / self.dz2)
                                  -  u * ((self.cl[i, 1:-1] - self.cl[i, :-2]) / self.dz) - \
                    3 * self.components[i].kf / (self.bc * self.rp) * (self.cl[i, 1:-1] - self.cp[i, 1:-1])) * dt
                    
            self.cl[i, -1] = self.cl[i, -1] - (u * (self.cl[i, -1] - self.cl[i, -2]) / (self.dz * 0.5)) * dt
                   
            self.c_out[i] = self.cl[i, -1]
         
        # Update particle concentrations
        for i in range(len(self.components)):
            self.cp[i][:] = self.cp[i][:] + \
                (3 * self.components[i].kf / (self.ep * self.rp * 1) * (self.cl[i][:] - self.cp[i][:])) * dt
        

class Component:
    def __init__ (self, kf=0):
        self.kf = kf
        
        
class Solution:
    def __init__ (self, components=[], concentrations=[]):
        self.components = components
        self.concentrations = concentrations
        

class Phase:
    def __init__ (self, inlet_solution, flowrate, t):
        self.inlet_solution = inlet_solution
        self.flowrate = flowrate
        self.t = t
    

class Experiment:
    def __init__(self, components, units, phases, dt):
        self.components = components
        self.units = units
        self.phases = phases
        self.t_elapsed = 0
        self.dt = dt
        
        self.current_phase = phases[0]
        
        self.t_solve = []
        self.c_solve = []
    
    def update_phases(self):
        for phase in self.phases:
            if self.t_elapsed <= phase.t:
                self.current_phase = phase
                return
            
    
    def step (self):
        self.update_phases()
        
        c0 = self.current_phase.inlet_solution.concentrations
        flowrate = self.current_phase.flowrate
        phase_length = self.current_phase.t
        


        c_previous = c0
        for unit in self.units:
            unit.step(c_previous, flowrate, self.dt)
            c_previous = unit.c_out
            
        self.c_solve.append(list(self.units[0].c_out))
        
        self.t_elapsed += self.dt
    
    
    def run (self):
        total_time = sum([i.t for i in self.phases])
        #self.units[0].components = self.components
        
        while self.t_elapsed <= total_time:
            self.step()
            self.t_solve.append(self.t_elapsed)
        
        for i in range(len(self.components)):
            plt.plot(self.t_solve, [x[i] for x in self.c_solve])
            
        plt.show()