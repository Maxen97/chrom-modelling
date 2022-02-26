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


#salt = Salt(charge=1)
salt = Component()
protein1 = Component()
protein2 = Component()

components = [salt, protein1, protein2]

buffer1 = Solution(components, [0.2, 1.8, 0.4])
buffer2 = Solution(components, [0.2, 0., 0.])

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
