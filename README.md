# Chromatography modelling

## About

`chrom-modelling` is a simple chromatography sandbox simulator.

Compared to other chromatography simulation tools, `chrom-modelling` will have built in support for parameter sensitivity analysis and an intuitive UI for easier adoption of multiple unit operations within the process.

Please be aware that development is ongoing.

## Usage

### Getting started
1. Download and install `numpy` and `matplotlib`.
2. Save the "chrom_modelling.py" script to a local repository.
3. Import the "chrom_modeling.py" script to your running script as seen in the example below.

### Example
The following example shows a model of a standard pulse injection experiment.

```python
import chrom_modelling as model # Import main script

# Define components (e.g., salts, biomolecules, additives)
salt = model.Component(
    kf=1e-5     # Film diffusion coefficient [m/s]
)

# Define solutions (e.g., load solutions, buffers)
load = model.Solution(
    components=[salt],
    concentrations=[0.0005] # [moles/m^3]
)
eluent = model.Solution(
    components=[salt],
    concentrations=[0.]     # [moles/m^3]
)

# Define model unit operations (e.g., columns, tubings, tanks)
column = model.LRMP(
    nz=50,      # Axial discretization [-]
    lz=0.2,     # Axial length [m]
    lr=0.01,    # Radial length [m]
    et=0.85,    # Total porosity [-]
    ec=0.35,    # Column porosity [-]
    rp=20e-6,   # Particle radius [m]
    dax=1e-5,   # Axial dispersion [m^2/s]
    components=[salt],
    initial_concentrations=[0.] # [moles/m^3]
)

# Define experiment phases (e.g., loading phase, elution phase)
injection_phase = model.Phase(
    inlet_solution=load,
    flowrate=1e-6, # [m^3/s]
    t=10 # [s]
)
elution_phase = model.Phase(
    inlet_solution=eluent,
    flowrate=1e-6, # [m^3/s]
    t=150 # [s]
)

# Create experiment
experiment = model.Experiment(
    components=[salt],
    units=[column],
    phases=[injection_phase, elution_phase],
    dt=0.01 # Seconds between time-steps [s]
)

experiment.run()
```