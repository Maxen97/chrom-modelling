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

buffer1 = Solution([(salt, 50)])
buffer2 = Solution([(salt, 1000)])
buffer3 = Solution([(salt, 50), (protein1, 200), (protein2, 160)])
"""


column = Unit(type="LRMP", length=0.10, Dax=0.00001, ax_disc=50)


init_sys_params = SystemParameters()
init_sys_params.inlet_concentration = 0.0
init_sys_params.u = 0.00001  #m3/s

apply = deepcopy(init_sys_params)
apply.inlet_concentration = 1.8

elute = deepcopy(init_sys_params)


experiment = Experiment()
experiment.system_parameters = init_sys_params
experiment.units = [column]
experiment.add_timesection(TimeSection(system_parameters=apply, time=6.1))
experiment.add_timesection(TimeSection(system_parameters=elute, time=5))
experiment.run()

