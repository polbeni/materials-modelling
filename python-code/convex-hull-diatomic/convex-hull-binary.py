# Pol Benítez Colominas, October 2025
# The University of Tokyo and Universitat Politècnica de Catalunya

# Determine the complex hull using pymatgen for a binary system

import matplotlib.pyplot as plt
import pandas as pd

from pymatgen.core import Composition, Element
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter


# Read the data from csv generated file
df = pd.read_csv("energies.csv")
#entries = [ComputedEntry(Composition(row["formula"]), row["energy_per_atom"]) for _, row in df.iterrows()]
entries = [
    ComputedEntry(
        Composition(row["formula"]),
        row["energy_per_atom"] * Composition(row["formula"]).num_atoms
    )
    for _, row in df.iterrows()
]

# Determine phase diagram using pymatgen tools
pd = PhaseDiagram(entries)

# Get the stable compounds, those with 0 energy above hull and save them in a file
stable = open('stable_compounds.txt', 'w')

for e in pd.stable_entries:
    print(f"Stable phase: {e.composition.reduced_formula}")
    stable.write(f'{e.composition.reduced_formula}\n')

stable.close()

# Get information for all the entries
data_x = []
data_y = []
data_z = []

all = open('all_compounds.txt', 'w')
all.write('formula    x_A    E_form_per_atom    E_above_hull\n')
for e in entries:
    comp = e.composition
    form_e = pd.get_form_energy_per_atom(e)          # formation energy per atom
    e_above = pd.get_e_above_hull(e)                 # energy above convex hull
    xA = comp.get_atomic_fraction("A")                   # fraction of A
    all.write(f'{comp.reduced_formula}    {xA}    {form_e}    {e_above}\n')

    data_x.append(xA)
    data_y.append(form_e)
    data_z.append(e_above)

all.close()

# Plot the results
plt.figure(figsize=(8, 6))
scatter = plt.scatter(data_x, data_y, c=data_z, cmap='viridis', s=100, edgecolor='k', vmin=0, vmax=.6)


cbar = plt.colorbar(scatter)
cbar.set_label('Z value')

plt.ylim(-.4 , 0.1)

plt.xlabel('X values')
plt.ylabel('Y values')
plt.title('Scatter plot colored by Z value')

plt.grid(True)
plt.show()
