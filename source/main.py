# -*- coding: utf-8 -*-
"""
Created on Sun Dec 26 17:59:00 2021

@author: Admin
"""

import numpy as np
import matplotlib.pyplot as plt
from copy import copy, deepcopy

from helper import *
from unit import Unit

"""
salt = Salt(charge=1)
protein1 = Protein()
protein2 = Protein()

components = [salt, protein1, protein2]

buffer1 = Solution(components, [50, 0, 0])
buffer2 = Solution(components, [1000, 0, 0])
buffer3 = Solution(components, [50, 200, 140])

phase1 = Phase(inlet_solution=buffer3,
               flowrate=0.00001,
               length=0.5)
phase2 = Phase(inlet_solution=buffer1,
               flowrate=0.00001,
               length=10)

experiment.add_phase(phase1)
experiment.add_phase(phase2)
"""




#salt = Salt(charge=1)
protein1 = Component()
protein2 = Component()

components = [protein1, protein2]

buffer1 = Solution(components, [1.8, 0.4])
buffer2 = Solution(components, [0., 0.])

phase1 = Phase(inlet_solution=buffer1,
               flowrate=0.00001,
               length=0.5)
phase2 = Phase(inlet_solution=buffer2,
               flowrate=0.00001,
               length=10)

column = Unit(type="LRMP", length=0.30, Dax=0.0001, ax_disc=50, components=components)

experiment = Experiment()
experiment.units = [column]
experiment.components=components
experiment.add_phase(phase1)
experiment.add_phase(phase2)
experiment.run()


"""
init_sys_params = SystemParameters()
init_sys_params.inlet_concentration = 0.0
init_sys_params.u = 0.00001  #m3/s

apply = deepcopy(init_sys_params)
apply.inlet_concentration = 1.8

elute = deepcopy(init_sys_params)


experiment = Experiment()
experiment.system_parameters = init_sys_params
experiment.units = [column]
experiment.add_timesection(TimeSection(system_parameters=apply, time=.1))
experiment.add_timesection(TimeSection(system_parameters=elute, time=8))
experiment.run()
"""
