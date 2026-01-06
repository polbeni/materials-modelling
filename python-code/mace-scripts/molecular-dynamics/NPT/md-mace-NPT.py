# Pol Benítez Colominas, Feb 2025 - Jan 2026
# The University of Tokyo and Universitat Politècnica de Catalunya

# Script to run a relaxation using ASE calculator and MACE machine-learning interatomic potentials
# Uses a NPT ensemble

import warnings

from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor

from ase.md.nptberendsen import NPTBerendsen
from ase.io.trajectory import Trajectory
from ase.io.vasp import write_vasp_xdatcar, write_vasp
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase import units

from mace.calculators import mace_mp
import torch
torch.serialization.add_safe_globals([slice])

#torch.backends.cuda.preferred_linalg_library("magma")

warnings.simplefilter('ignore')

# Choose the device you want to use cpu or cuda
device = 'cuda' #'cuda'

# Chose the model (default: 'large')
model = 'mace-omat-0-medium.model'

# Read the system (with VASP structure format)
crystal_structure = Structure.from_file('SPOSCAR')
ase_adaptor = AseAtomsAdaptor()
atoms = ase_adaptor.get_atoms(crystal_structure)

# Load the MACE calculator
atoms.calc = mace_mp(model=model, device=device)

# Set up the molecular dynamics parameters
temperature = 400      # in K
timestep = 1           # in fs
n_steps = 2000        # number of steps in the simulation
nblock = 10            # interval of saved steps
friction = 5e-3        # in fs^-1, a friction coefficient

# To ensure we are using the correct units in time dependent variables
# This is important since time units are given in other units, read: https://wiki.fysik.dtu.dk/ase/ase/units.html
timestep = timestep * units.fs

# Initialize the velocities using a Maxwell-Boltzmann distribution
# It is also possible to use phononic data to initialize the velocitities (read: https://wiki.fysik.dtu.dk/ase/ase/md.html#velocity-distributions)
MaxwellBoltzmannDistribution(atoms, temperature_K=temperature)

# Perform the molecular dynamics (here using Langevin, since we want NVT ensemble, 
# other method can be found in ASE documentation: https://wiki.fysik.dtu.dk/ase/ase/md.html)
dyn = NPTBerendsen(
    atoms=atoms,                              # the list of atoms
    timestep=timestep,                        # the time step in ASE time units
    temperature_K=temperature,                # the desired temperature, in Kelvin
    pressure_au = 0.0,                        # the desired pressure, in atomic units (eV/Å^3), in this case no pressure
    taut=100 * units.fs,                      # time constant for Berendsen temperature coupling in ASE time units (large value for less agressive MD)
    taup=1000 * units.fs,                     # time constant for Berendsen pressure coupling (large value for less agressive MD)
    compressibility_au=5.5e-7 / units.bar,    # the compressibility of the material, in atomic units (Å^3/eV)
    trajectory='run.traj',                    # attach trajectory object
    logfile='run.log',                        # file with log data
    loginterval=nblock                        # interval of saved steps
)

dyn.run(n_steps)

# Save the relaxed structure (with VASP structure format)
bec =  Trajectory('run.traj')
bec_new_cont = write_vasp('CONTCAR', atoms=atoms, direct=True)
bec_new_xdat = write_vasp_xdatcar('XDATCAR', bec, label=None)
