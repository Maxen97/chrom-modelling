import matplotlib.pyplot as plt
from chrom_modelling import model

# Define components
salt = Component(
    type = "salt",

)
protein1 = Component()
protein2 = Component()

components = [salt, protein1, protein2]


# Define solutions
buffer1 = Solution(components, [0.2, 1.8, 0.4])
buffer2 = Solution(components, [0.2, 0., 0.])


# Define system
column = Unit(type="LRMP", length=0.30, Dax=0.0001, ax_disc=50, components=components)
column.set_init_concentrations([0.2, 0.0, 0.0])


# Add phases
phase1 = Phase(inlet_solution=buffer1,
               flowrate=0.00001,
               length=0.5)
phase2 = Phase(inlet_solution=buffer2,
               flowrate=0.00001,
               length=10)


# Create experiment
experiment = Experiment()
experiment.units = [column]
experiment.components=components
experiment.add_phase(phase1)
experiment.add_phase(phase2)
experiment.run()