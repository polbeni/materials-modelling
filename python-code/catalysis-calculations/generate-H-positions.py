# Pol Benítez Colominas, July 2024
# Universitat Politècnica de Catalunya

# This script generates POSCAR files with hydrogen atoms in different positions
# of a slab, in order to perform catalysis calculations

import os
import shutil
from copy import deepcopy

from pymatgen.core import Structure, Element

# Open the slab structure
structure = Structure.from_file('POSCAR')

# Define the atom (or atoms) to add
new_atom = 'H'

# Define the distance of the atom to the slab
distance_atom = 1.2 # in angstrom
z_axis_length = structure.lattice.c
distance_frac = 1.2/z_axis_length

# Create a list of positions for the atom
# because of symmetry we only need to consider few positions (manually)
positions = []
# top slab (there is no top-bottom slab symmetry)
positions.append([0.5, 0.5, 0.73159 + distance_frac]) # above bromine
positions.append([0.25, 0.5, 0.73159 + distance_frac]) # above middle silver
positions.append([0.25, 0.25, 0.73159 + distance_frac]) # above cornen silver
positions.append([0.375, 0.5, 0.73159 + distance_frac]) 
positions.append([0.375, 0.375, 0.73159 + distance_frac])
# bottom slab
positions.append([0.5, 0.5, 0.26841 - distance_frac]) # above bromine
positions.append([0.25, 0.5, 0.26841 - distance_frac]) # above middle silver
positions.append([0.25, 0.25, 0.26841 - distance_frac]) # above cornen silver
positions.append([0.375, 0.5, 0.26841 - distance_frac]) 
positions.append([0.375, 0.375, 0.26841 - distance_frac])

# Generate the structures
if os.path.exists('slab_hydrogen_structures'):
    shutil.rmtree('slab_hydrogen_structures')
os.mkdir('slab_hydrogen_structures')

iteration = 1
path_save = 'slab_hydrogen_structures/POSCAR-'
for position in positions:
    slab = deepcopy(structure)

    slab.append(new_atom, position)

    slab.to(fmt='poscar', filename=path_save+str(iteration).zfill(3))

    iteration = iteration + 1