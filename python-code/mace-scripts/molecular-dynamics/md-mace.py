# Pol Benítez Colominas, February 2025
# Universitat Politècnica de Catalunya

# Script to run a relaxation using ASE calculator and MACE machine-learning interatomic potentials
# The returned values are the system energy, the atoms forces and the system stress tensor

import warnings

from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor

from ase.md import Langevin
from ase.io.trajectory import Trajectory
from ase.io.vasp import write_vasp_xdatcar, write_vasp
from ase import units

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
n_steps = 100          # number of steps in the simulation
nblock = 10            # interval of saved steps
friction = 5e-3        # in fs^-1, a friction coefficient

# To ensure we are using the correct units in time dependent variables
# This is important since time units are given in other units, read: https://wiki.fysik.dtu.dk/ase/ase/units.html
timestep = timestep * units.fs
friction = friction / units.fs
# IMPORTANT: for temperature we do not need to do this since we are already using tempereature_K when 
# running the md in Langevin object

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


########### MD ALGORITHMS ###########

# NVE ensemble: VelocityVerlet
# Read the documentation: https://wiki.fysik.dtu.dk/ase/ase/md.html#constant-nve-simulations-the-microcanonical-ensemble

# NVT ensemble: Langevin, NoseHooverChainNVT, Bussi (the following algorithms are not recommended: Andersen, NVTBerendsen)
# Read the documentation: https://wiki.fysik.dtu.dk/ase/ase/md.html#constant-nvt-simulations-the-canonical-ensemble

# NPT ensemble (ASE lacks a good NPT algorithm): NPTBerendsen, NPT
# Read the documentation: https://wiki.fysik.dtu.dk/ase/ase/md.html#constant-npt-simulations-the-isothermal-isobaric-ensemble

#####################################