# Pol Benítez Colominas, February 2025
# Universitat Politècnica de Catalunya

# Script to run a relaxation using ASE calculator and MACE machine-learning interatomic potentials
# The returned values are the system energy, the atoms forces and the system stress tensor

import warnings

from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor

from ase.optimize import BFGS, FIRE
from ase.io.trajectory           import Trajectory
from ase.io.vasp                 import write_vasp
from ase.constraints import ExpCellFilter

from mace.calculators import mace_mp
import torch
torch.serialization.add_safe_globals([slice])

warnings.simplefilter('ignore')

# Read the system (with VASP structure format)
crystal_structure = Structure.from_file('POSCAR')
ase_adaptor = AseAtomsAdaptor()
atoms = ase_adaptor.get_atoms(crystal_structure)

# Load the MACE calculator
atoms.calc = mace_mp(model='large', device='cpu')

# Allow lattice parameters to change
atoms_filter = ExpCellFilter(atoms)

# Set up the relaxation parameters
fmax = 0.05 # max force allowed
smax = 200 # max number of steps for the relaxation

# Perform the relaxation (here using BFGS other method can be found in 
# ASE documentation: https://wiki.fysik.dtu.dk/ase/ase/optimize.html)
dyn = BFGS(atoms_filter, trajectory='run.traj', logfile='run.log')
dyn.run(fmax=fmax, steps=smax)

# Save the relaxed structure (with VASP structure format)
bec =  Trajectory('run.traj')
bec_new_cont = write_vasp('CONTCAR', atoms=atoms, direct=True)