# Pol Benítez Colominas, November 2023
# Universitat Politècnica de Catalunya

# Code to compute phonons in harmonic approximation using hiPhive results

import numpy as np
import matplotlib.pyplot as plt

from ase.io.vasp import read_vasp

from ase import Atoms
from ase.build import bulk
from hiphive import ForceConstantPotential

from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms

# define the parameters
dim = 2 # supercell size
Nq = 51  # number of q-points along each segment of the path through the BZ

# create phonopy object for the desired material
prim = read_vasp(file='POSCAR')
atoms_phonopy = PhonopyAtoms(symbols=prim.get_chemical_symbols(),
                             scaled_positions=prim.get_scaled_positions(),
                             cell=prim.cell)
phonopy = Phonopy(atoms_phonopy, supercell_matrix=dim*np.eye(3),
                  primitive_matrix=None)

# read the force constant potential obtained with hiPhive
fcp = ForceConstantPotential.read('Ag3SI.fcp')

# find the force constant
supercell = phonopy.get_supercell()
supercell = Atoms(cell=supercell.cell, numbers=supercell.numbers, pbc=True,
                  scaled_positions=supercell.get_scaled_positions())
fcs = fcp.get_force_constants(supercell)

phonopy.set_force_constants(fcs.get_fc_array(order=2))

# determine the bands for the desired k-path
def get_band(q_start, q_stop, N):
    """ Return path between q_start and q_stop """
    return np.array([q_start + (q_stop-q_start)*i/(N-1) for i in range(N)])


G2X = get_band(np.array([0, 0, 0]), np.array([0, 0.5, 0]), Nq)
X2M = get_band(np.array([0, 0.5, 0]), np.array([0.5, 0.5, 0]), Nq)
M2G = get_band(np.array([0.5, 0.5, 0]), np.array([0, 0, 0]), Nq)
G2R = get_band(np.array([0, 0, 0]), np.array([0.5, 0.5, 0.5]), Nq)
R2X = get_band(np.array([0.5, 0.5, 0.5]), np.array([0, 0.5, 0]), Nq)
R2M = get_band(np.array([0.5, 0.5, 0.5]), np.array([0.5, 0.5, 0]), Nq)

bands = [G2X, X2M, M2G, G2R, R2X, R2M]

phonopy.set_band_structure(bands)
qvecs, qnorms, freqs, _ = phonopy.get_band_structure()

print(freqs)