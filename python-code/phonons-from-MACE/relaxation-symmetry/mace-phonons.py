# Pol Benítez Colominas, February 2025
# Universitat Politècnica de Catalunya

# It computes the phonon dispersion (anharmonic) using phonopy from the forces determined with MACE
# The workingflow of the script is:
#       1. Relax the structure with MACE
#       2. Generate the distorted structures with phonopy for the desired supercell
#       3. Compute the forces with MACE for all the distorted structures and generate FORCE_SETS
#       4. Calculate the phonon dispersion (harmonic) for the desired path using phonopy
#       5. Perform a molecular dynamics at a given temperature with MACE
#       6. Renormalize phonon frequencies with temperature using DynaPhoPy
#       7. Calculate the phonon dispersion (anharmonic) using phonopy

# In order to run the script we need to provide the following input files: POSCAR, phonopy.conf
# We need to provide the following variables: poscar_name, dimension_supercell, device, model, temperature
#       poscar_name -> name of the initial POSCAR structure (not relaxed)
#       dimension_supercell -> dimension of the supercell in str format: "2 2 2"
#       device -> cpu or cuda for mace calculations
#       model -> MACE model to use (default: large)
#       temperature -> temperature of the molecular dynamics simulation to renormalize phonon frequencies


import subprocess
import yaml
import warnings

from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor

from ase.optimize import BFGS, FIRE
from ase.md import Langevin
from ase.io.trajectory import Trajectory
from ase.io.vasp import write_vasp, write_vasp_xdatcar
from ase.constraints import ExpCellFilter, FixAtoms
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase import units

from mace.calculators import mace_mp
import torch
torch.serialization.add_safe_globals([slice])

warnings.simplefilter('ignore')


### Variables definition
poscar_name = 'original_POSCAR'
dimension_supercell = '2 2 2'
device = 'cuda'
model = 'large'
temperature = 200


### Read the original structure, relax it and save it with the name POSCAR
crystal_structure = Structure.from_file(poscar_name) # read the initial structure
ase_adaptor = AseAtomsAdaptor()
atoms = ase_adaptor.get_atoms(crystal_structure)

atoms.calc = mace_mp(model=model, device=device) # load the MACE calculator

constraint = FixAtoms(indices=[atom.index for atom in atoms]) # fix the atoms positions in the relaxation
atoms.set_constraint(constraint)

mask = [1, 1, 1, 0, 0, 0]
atoms_filter = ExpCellFilter(atoms, mask=mask) # allow lattice parameters to change

dyn = BFGS(atoms_filter) # relax the structure
dyn.run(fmax=0.001, steps=200)

write_vasp('POSCAR', atoms=atoms, direct=True) # save the relaxed structure


### Use phonopy to generate the distorted structures
command = f'phonopy -d --dim="{dimension_supercell}"'  
result = subprocess.run(command, shell=True, capture_output=True, text=True)

### Save the number of atoms in the supercell and the number of distorted structures
with open('phonopy_disp.yaml', 'r') as f:
    data = yaml.safe_load(f)

atom_info = data['primitive_cell']['points']
num_atoms_cell = len(atom_info)
num_atoms = num_atoms_cell # in the supercell
for component in range(len(dimension_supercell.split())): # diagonal supercell assumed
    num_atoms = num_atoms*int(dimension_supercell.split()[component])

num_distorted_struc = len(data['displacements'])


### Compute the forces with MACE for all the distorted structures
forces_array = []
for dist in range(num_distorted_struc):
    crystal_structure = Structure.from_file(f'POSCAR-{str(dist + 1).zfill(3)}') # read the distorted structure
    ase_adaptor = AseAtomsAdaptor()
    atoms = ase_adaptor.get_atoms(crystal_structure)

    atoms.calc = mace_mp(model=model, device=device) # load mace calculator 

    forces = atoms.get_forces() # get the forces
    forces_array.append(forces)


### Generate the FORCE_SETS file
force_sets = open('FORCE_SETS', 'w')

force_sets.write(f'{num_atoms}\n')
force_sets.write(f'{num_distorted_struc}\n')
force_sets.write('\n')

for dist in range(num_distorted_struc):
    atom_ind = data['displacements'][dist]['atom']
    force_vec = data['displacements'][dist]['displacement']

    force_sets.write(f'{atom_ind}\n')
    force_sets.write(f'  {force_vec[0]}   {force_vec[1]}   {force_vec[2]}\n')
    for atom in range(num_atoms):
        force_sets.write(f'     {forces_array[dist][atom][0]}   {forces_array[dist][atom][1]}   {forces_array[dist][atom][2]}\n')
    force_sets.write('\n')

force_sets.close()


### Calculate phonon dispersions (harmonic) with phonopy
command = 'phonopy -ps phonopy-harmonic.conf'  
result = subprocess.run(command, shell=True, capture_output=True, text=True)

command = 'mv band_dos.pdf band_dos-harmonic.pdf'  
result = subprocess.run(command, shell=True, capture_output=True, text=True)
command = 'mv band.yaml band-harmonic.yaml'  
result = subprocess.run(command, shell=True, capture_output=True, text=True)
command = 'mv mesh.yaml mesh-harmonic.yaml'  
result = subprocess.run(command, shell=True, capture_output=True, text=True)
command = 'mv total_dos.dat total_dos-harmonic.dat'  
result = subprocess.run(command, shell=True, capture_output=True, text=True)


### Run a molecular dynamics simulation at the desired temperature
crystal_structure = Structure.from_file('SPOSCAR') # read the supercell structure
ase_adaptor = AseAtomsAdaptor()
atoms = ase_adaptor.get_atoms(crystal_structure)

atoms.calc = mace_mp(model=model, device=device) # load the MACE calculator

timestep = 1             # in fs
n_steps = 40000          # number of steps in the simulation
nblock = 10              # interval of saved steps
friction = 5e-3          # in fs^-1, a friction coefficien

timestep = timestep * units.fs # to ensure we are using the correct units in time dependent variables
friction = friction / units.fs

MaxwellBoltzmannDistribution(atoms, temperature_K=temperature) # initialize the velocities

dyn = Langevin(
    atoms=atoms,                 # the list of atoms
    timestep=timestep,           # the time step in ASE time units
    temperature_K=temperature,   # the desired temperature, in Kelvin
    friction=friction,           # a friction coefficient in inverse ASE time units
    trajectory='run.traj',       # attach trajectory object
    logfile='run.log',           # file with log data
    loginterval=nblock           # interval of saved steps
)

dyn.run(n_steps) # run the md

bec =  Trajectory('run.traj') # save the results in XDATCAR file
bec_new_xdat = write_vasp_xdatcar('XDATCAR', bec, label=None)


### Renormalize the phonon frequencies with temperature using dynaphopy
time_step = nblock * timestep * 1e3 # in ps
command = 'dynaphopy input XDATCAR -ts 0.01 -sfc FORCE_CONSTANTS -psm 2'  
result = subprocess.run(command, shell=True, capture_output=True, text=True)


### Calculate phonon dispersions (anharmonic) with phonopy
command = 'phonopy -ps phonopy-anharmonic.conf'  
result = subprocess.run(command, shell=True, capture_output=True, text=True)

command = 'mv band_dos.pdf band_dos-anharmonic.pdf'  
result = subprocess.run(command, shell=True, capture_output=True, text=True)
command = 'mv band.yaml band-anharmonic.yaml'  
result = subprocess.run(command, shell=True, capture_output=True, text=True)
command = 'mv mesh.yaml mesh-anharmonic.yaml'  
result = subprocess.run(command, shell=True, capture_output=True, text=True)
command = 'mv total_dos.dat total_dos-anharmonic.dat'  
result = subprocess.run(command, shell=True, capture_output=True, text=True)
