import matplotlib.pyplot as plt
from chrom_modelling import model

# Define components (e.g., salts, biomolecules, additives)
salt = model.Component(
    type = "salt",

)

# Define solutions (e.g., load solutions, buffers)
load = model.Solution(
    components=[salt],
    concentrations=[0.4] # [M]
)
eluent = model.Solution(
    components=[salt],
    concentrations=[0]
)

# Define model units (e.g., columns, tubings, tanks)
column = model.Unit(
    type="LRMP", # e.g., Lumped-Rate Model with Pores
    nz=50, # Axial discretization [-]
    lz=0.2, # Axial length [m]
    lr=0.05, # Radial length [m]
    components=[salt],
    init_concentrations=[0.4] # [M]
)

# Define experiment phases (e.g., loading phase, elution phase)
injection_phase = model.Phase(
    inlet_solution=load,
    flowrate=0.001, # [L/s]
    v_end=0.01 # [L]
)
elution_phase = model.Phase(
    inlet_solution=eluent,
    flowrate=0.001,
    v_end=0.1
)


# Create experiment
experiment = model.Experiment(
    components=[salt],
    units=[column],
    phases=[injection_phase, elution_phase]
)

experiment.run()