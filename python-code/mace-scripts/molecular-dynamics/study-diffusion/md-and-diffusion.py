# Pol Benítez Colominas, February 2025
# Universitat Politècnica de Catalunya

# Script to run a relaxation using ASE calculator and MACE machine-learning interatomic potentials
# and obtain the diffusion coefficients for all the atom in the system (using the Einstein equation)

import warnings

from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor

from ase.md import Langevin
from ase.io.trajectory import Trajectory
from ase.io.vasp import write_vasp_xdatcar, write_vasp
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase import units
from ase.md.analysis import DiffusionCoefficient

from mace.calculators import mace_mp
import torch
torch.serialization.add_safe_globals([slice])

warnings.simplefilter('ignore')

# Choose the device you want to use cpu or cuda
device = 'cpu' #'cuda'

# Chose the model (default: 'large')
model = 'large'

# Read the system (with VASP structure format)
crystal_structure = Structure.from_file('SPOSCAR')
ase_adaptor = AseAtomsAdaptor()
atoms = ase_adaptor.get_atoms(crystal_structure)

# Load the MACE calculator
atoms.calc = mace_mp(model=model, device=device)

# Set up the molecular dynamics parameters
temperature = 300      # in K
timestep = 1           # in fs
n_steps = 10          # number of steps in the simulation
nblock = 5            # interval of saved steps
friction = 5e-3        # in fs^-1, a friction coefficient

# To ensure we are using the correct units in time dependent variables
# This is important since time units are given in other units, read: https://wiki.fysik.dtu.dk/ase/ase/units.html
timestep = timestep * units.fs
friction = friction / units.fs
# IMPORTANT: for temperature we do not need to do this since we are already using tempereature_K when 
# running the md in Langevin object

# Initialize the velocities using a Maxwell-Boltzmann distribution
# It is also possible to use phononic data to initialize the velocitities (read: https://wiki.fysik.dtu.dk/ase/ase/md.html#velocity-distributions)
MaxwellBoltzmannDistribution(atoms, temperature_K=temperature)

# Perform the molecular dynamics (here using Langevin, since we want NVT ensemble, 
# other method can be found in ASE documentation: https://wiki.fysik.dtu.dk/ase/ase/md.html)
dyn = Langevin(
    atoms=atoms,                 # the list of atoms
    timestep=timestep,           # the time step in ASE time units
    temperature_K=temperature,   # the desired temperature, in Kelvin
    friction=friction,           # a friction coefficient in inverse ASE time units
    trajectory='run.traj',       # attach trajectory object
    logfile='run.log',           # file with log data
    loginterval=nblock           # interval of saved steps
)

dyn.run(n_steps)

# Save the relaxed structure (with VASP structure format)
bec =  Trajectory('run.traj')
bec_new_cont = write_vasp('CONTCAR', atoms=atoms, direct=True)
bec_new_xdat = write_vasp_xdatcar('XDATCAR', bec, label=None)

# Compute the diffusion coefficients from the trajectories
diffusion_coefficient = DiffusionCoefficient(
    traj=bec,                             # trajectory object
    timestep=timestep * nblock,           # timestep between saved frames (in ASE units)
    atom_indices=list(range(len(atoms)))  # indices of atoms to include
)

diffusion_coefficient.calculate(
    ignore_n_images=0,        # number of initial images to ignore (e.g., equilibration)
    number_of_segments=1      # number of segments for statistical analysis
)

slopes, std = diffusion_coefficient.get_diffusion_coefficients()

# Save the Diffusion coefficient for each specie in a output file
output_file = open('diffusion.output', 'w')
output_file.write('Specie      D [Ang^2/fs]      std dev [Ang^2/fs]\n')
for i, (slope, stdev) in enumerate(zip(slopes, std)):
    output_file.write(f'{diffusion_coefficient.types_of_atoms[i]}      {slope}      {stdev}\n')
output_file.close()