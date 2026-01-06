# Pol Benítez Colominas, January 2026
# The University of Tokyo and Universitat Politècnica de Catalunya

# Read a run.traj file an extract the lattice change


import numpy as np
import pandas as pd

from ase.io import Trajectory


# Load the trajectory
traj = Trajectory('run.traj')

# Save the lattice values and the volume of the box at each step
steps = []
a_vals, b_vals, c_vals = [], [], []
volumes = []

for i, atoms in enumerate(traj):
    steps.append(i)
    lengths = atoms.cell.lengths()
    a_vals.append(lengths[0])
    b_vals.append(lengths[1])
    c_vals.append(lengths[2])
    volumes.append(atoms.get_volume())

# Save the results in a csv file
data = {'step': steps, 'a': a_vals, 'b': b_vals, 'c': c_vals, 'volume': volumes}
df = pd.DataFrame(data)
df.to_csv('lattice_evolution.csv', index=False)
