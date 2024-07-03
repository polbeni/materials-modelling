# Pol Benítez Colominas, date of creation: 2023/09/07, date of last modification: 2023/09/07

# Code to run molecular dynamics simulations with m3gnet model.
# Provide a POSCAR file and input file with the different variables for
# the simulation, and this code will generate a XDATCAR with the results.

import os

from pymatgen.core import Structure, Lattice
from m3gnet.models import MolecularDynamics

import ase
from ase.io.trajectory import Trajectory
from ase.io import read, write
from ase.io.vasp import write_vasp_xdatcar

# Read the POSCAR file
crystal_structure = Structure.from_file("POSCAR")

# Default values if no inputs file provided
tbeg = 1000
mdalgo = 'nvt'
potim = 1
nblock = 100
nsw = 1000

# Read the parameters from inputs file
if os.path.exists('inputs'):
    inputs = open('inputs', "r")

    variables = [0, 0, 0, 0, 0]

    for x in range(len(variables)):
        line = inputs.readline()

        variables[x] = line.split()[2]

    inputs.close()

    tbeg = float(variables[0])
    mdalgo = variables[1]
    potim = float(variables[2])
    nblock = float(variables[3])
    nsw = float(variables[4])

# Create the simulation model
md = MolecularDynamics(
    atoms=crystal_structure, #crsytal structure
    temperature=tbeg,  # temperature in K
    ensemble=mdalgo,  # type of ensemble
    timestep=potim, # time for each step in fs
    trajectory="mo.traj",  # save trajectory to mo.traj
    logfile="mo.log",  # log file for MD
    loginterval=nblock,  # interval for record the log
)

# Run the model
md.run(steps=nsw) 

# Save the results in a XDATCAR type file
bec =  Trajectory("mo.traj")
bec_new_xdat = write_vasp_xdatcar("XDATCAR", bec, label=None)