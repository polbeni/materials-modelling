# Pol Benítez Colominas, February 2025
# Universitat Politècnica de Catalunya

# It computes the phonon dispersion (harmonic) using phonopy from the forces determined with MACE
# The workingflow of the script is:
#       1. Relax the structure with MACE
#       2. Generate the distorted structures with phonopy for the desired supercell
#       3. Compute the forces with MACE for all the distorted structures and generate FORCE_SETS
#       4. Calculate the phonon dispersion for the desired path using phonopy

# In order to run the script we need to provide the following input files: POSCAR, phonopy.conf
# We need to provide the following variables: poscar_name, dimension_supercell, device, model
#       poscar_name -> name of the initial POSCAR structure (not relaxed)
#       dimension_supercell -> dimension of the supercell in str format: "2 2 2"
#       device -> cpu or cuda for mace calculations
#       model -> MACE model to use (default: large)

import subprocess
import yaml
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


### Variables definition
poscar_name = 'original_POSCAR'
dimension_supercell = '2 2 2'
device = 'cpu'
model = 'large'


### Read the original structure, relax it and save it with the name POSCAR
crystal_structure = Structure.from_file(poscar_name) # read the initial structure
ase_adaptor = AseAtomsAdaptor()
atoms = ase_adaptor.get_atoms(crystal_structure)

atoms.calc = mace_mp(model=model, device=device) # load the MACE calculator

atoms_filter = ExpCellFilter(atoms) # allow lattice parameters to change

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

    atoms.calc = mace_mp(model=model, device=device) # calculate 

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


### Calculate phonon dispersions with phonopy
command = 'phonopy -ps phonopy.conf'  
result = subprocess.run(command, shell=True, capture_output=True, text=True)
