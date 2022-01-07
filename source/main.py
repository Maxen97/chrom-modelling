# -*- coding: utf-8 -*-
"""
Created on Sun Dec 26 17:59:00 2021

@author: Admin
"""

import numpy as np
import matplotlib.pyplot as plt
from helper import *
from copy import copy, deepcopy

inlet = Unit(type="INLET")
column = Unit(type="COLUMN", length=3, Dax=0, ax_disc = 50)
outlet = Unit(type="OUTLET")


init_sys_params = SystemParameters()
init_sys_params.inlet_concentration = 0.0
init_sys_params.u = 4.0

apply = deepcopy(init_sys_params)
apply.inlet_concentration = 1.8

elute = deepcopy(init_sys_params)


experiment = Experiment()
experiment.system_parameters = init_sys_params
experiment.units = [inlet, column, outlet]
experiment.add_timesection(TimeSection(system_parameters=apply, time=0.1))
experiment.add_timesection(TimeSection(system_parameters=elute, time=1))
experiment.run()

