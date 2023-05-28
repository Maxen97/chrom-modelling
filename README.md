# Chromatography modeling

## About

`chrom-modelling` is a simple chromatography sandbox simulator.

Please be aware that development is ongoing.

## Usage

### Getting started
1. Copy the "models" folder a local repository.
2. Import the "models" folder to your running script as seen in the example below.

### Example
The following example shows a model of a standard pulse injection experiment.

```python
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
    concentrations=[0] # [M]
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
    t=10 # [s]
)
elution_phase = model.Phase(
    inlet_solution=eluent,
    flowrate=0.001, # [L/s]
    t=500 # [s]
)


# Create experiment
experiment = model.Experiment(
    components=[salt],
    units=[column],
    phases=[injection_phase, elution_phase]
)

result = model.solve(
    experiments=[experiment],
    dt=0.1 # [s]
)

fig = plt.plot(result.x, result.y)
fig.show()
```