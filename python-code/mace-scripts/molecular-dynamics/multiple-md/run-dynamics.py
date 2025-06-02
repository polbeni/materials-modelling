# Pol Benítez Colominas, June 2025
# Universitat Politècnica de Catalunya

# Run dynamics for different solid solutions and temperatures

import os
import shutil

import warnings

from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor

from ase.md import Langevin
from ase.io.trajectory import Trajectory
from ase.io.vasp import write_vasp_xdatcar, write_vasp
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase import units

from mace.calculators import mace_mp
import torch
torch.serialization.add_safe_globals([slice])

warnings.simplefilter('ignore')

# Set up the molecular dynamics parameters
timestep = 1           # in fs
n_steps = 60000        # number of steps in the simulation
nblock = 50            # interval of saved steps
friction = 5e-3        # in fs^-1, a friction coefficient
# Use the correct units
timestep = timestep * units.fs
friction = friction / units.fs

# Choose the device you want to use cpu or cuda
device = 'cuda' #'cpu'

# Chose the model (default: 'large')
model = 'mace-mpa-0-medium.model' #'large'

# Define the range of temperatures and materials
temperatures = [100, 200, 300, 400, 500, 600] # in K
ss_names_list = ['Ag3SBr', 'Ag3SBr0_125I0_875', 'Ag3SBr0_25I0_75', 'Ag3SBr0_375I0_625', 'Ag3SBr0_5I0_5', 'Ag3SBr0_625I0_375', 
                 'Ag3SBr0_75I0_25', 'Ag3SBr0_875I0_125', 'Ag3SI']

# Generate a path for the results
if os.path.exists('results'):
    shutil.rmtree('results')
os.mkdir('results')

# Perform the calculations for each solid solution and temperature
for material in ss_names_list:
    path_results = 'results/' + material
    if os.path.exists(path_results):
        shutil.rmtree(path_results)
    os.mkdir(path_results)

    for temp in temperatures:
        path_calculation = path_results + '/' + str(temp)
        if os.path.exists(path_calculation):
            shutil.rmtree(path_calculation)
        os.mkdir(path_calculation)

        # Load the structure
        path_struc = 'structures/supercells/' + material + '.vasp'
        crystal_structure = Structure.from_file(path_struc)
        ase_adaptor = AseAtomsAdaptor()
        atoms = ase_adaptor.get_atoms(crystal_structure)

        # Load the MACE calculator
        atoms.calc = mace_mp(model=model, device=device)

        # Define the value of the temperature
        temperature = temp

        # Initialize the velocities using a Maxwell-Boltzmann distribution
        MaxwellBoltzmannDistribution(atoms, temperature_K=temperature)

        # Perform the molecular dynamics
        dyn = Langevin(
            atoms=atoms,                                     # the list of atoms
            timestep=timestep,                               # the time step in ASE time units
            temperature_K=temperature,                       # the desired temperature, in Kelvin
            friction=friction,                               # a friction coefficient in inverse ASE time units
            trajectory=path_calculation + '/run.traj',       # attach trajectory object
            logfile=path_calculation + '/run.log',           # file with log data
            loginterval=nblock                               # interval of saved steps
        )

        dyn.run(n_steps)

        # Save the relaxed structure (with VASP structure format)
        bec =  Trajectory(path_calculation + '/run.traj')
        bec_new_cont = write_vasp(path_calculation + '/CONTCAR', atoms=atoms, direct=True)
        bec_new_xdat = write_vasp_xdatcar(path_calculation + '/XDATCAR', bec, label=None)
