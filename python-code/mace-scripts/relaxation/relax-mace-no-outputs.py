# Pol Benítez Colominas, March 2025
# Universitat Politècnica de Catalunya

# Script to run a relaxation using ASE calculator and MACE machine-learning interatomic potentials
# It does not use output file (run.traj)

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

# Choose the device you want to use cpu or cuda
device = 'cpu' #'cuda'

# Chose the model (default: 'large')
model = 'large'

# Read the system (with VASP structure format)
crystal_structure = Structure.from_file('POSCAR')
ase_adaptor = AseAtomsAdaptor()
atoms = ase_adaptor.get_atoms(crystal_structure)

# Load the MACE calculator
atoms.calc = mace_mp(model=model, device=device)

# Allow lattice parameters to change
atoms_filter = ExpCellFilter(atoms)

# Set up the relaxation parameters
fmax = 0.05 # max force allowed
smax = 200 # max number of steps for the relaxation

# Perform the relaxation (here using BFGS other method can be found in 
# ASE documentation: https://wiki.fysik.dtu.dk/ase/ase/optimize.html)
dyn = BFGS(atoms_filter, logfile='run.log')
dyn.run(fmax=fmax, steps=smax)

# Get the final relaxed structure
relaxed_atoms = atoms_filter.atoms
relaxed_structure = AseAtomsAdaptor().get_structure(relaxed_atoms)

# Save in VASP format
relaxed_structure.to(filename="CONTCAR", fmt="poscar")

# Get and print the final energy of the structure
final_energy = relaxed_atoms.get_potential_energy()

num_atoms = len(relaxed_atoms)
energy_per_atom = final_energy / num_atoms

print(f"Final energy per atom: {energy_per_atom:.4f} eV/atom")
