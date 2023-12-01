# Pol Benítez Colominas, November 2023
# Universitat Politècnica de Catalunya

# Code to generate the force constants with hiPhive from the rattled
# structures results from VASP

import os
import sys

from ase.io import read, write
from hiphive import ClusterSpace, StructureContainer, ForceConstantPotential
from hiphive.utilities import prepare_structures
from trainstation import Optimizer

# parameters
cutoffs = [4.5, 4.0, 4.0]               # cutoffs for pairs, triplets and quadruplets (in angstrom)
directory_path = 'data/Ag3SI-Pm-3m'     # path of the VASP results
optimization_method = 'least-squares'   # name of the desired optimization method employed

# read the structure input files
prim = read('prim.extxyz') # primitive cell
atoms_ideal = read('supercell_ideal.extxyz') # ideal supercell

# create the cluster space for the defined cutoffs
cs = ClusterSpace(prim, cutoffs)

os.mkdir('log_files')

with open('log_files/1-cluster-space', 'w') as f: # save the information in a file
    print('', cs, file=f)

original_stdout = sys.stdout
with open('log_files/2-cluster-space-orbits', 'w') as file: # save the information in a file
    sys.stdout = file
    cs.print_orbits()
sys.stdout = original_stdout

# define a function to get the different directories with rattled structures results
def get_directories_sorted(path, starts_with):
    directories = [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d)) and d.startswith(starts_with)]
    
    def extract_number(directory):
        return int(directory.split('-')[-1]) 

    sorted_directories = sorted(directories, key=extract_number)
    
    return sorted_directories

# generate a list with the name of these directories
prefix = 'disp'

disp_direct = get_directories_sorted(directory_path, prefix)

# generate the rattled file, i.e., a file that contains energy, positions and
# forces for each of the rattled structures
rattled_file = open('rattled_structures.extxyz', "w")

for disp in disp_direct:
    path_structure = directory_path + '/' + disp + '/vasprun.xml'
    rattled = read(path_structure, format='vasp-xml')

    write('actual_rattled.extxyz', rattled)

    #rattled_actual = open('actual_rattled.extxyz', "r")
    with open('actual_rattled.extxyz', 'r') as file:
        for line in file:
            rattled_file.write(line)

os.remove('actual_rattled.extxyz')
rattled_file.close()

# read the rattled strucutres file
rattled_structures = read('rattled_structures.extxyz', index=':')

# generate the structure container
structures = prepare_structures(rattled_structures, atoms_ideal)
sc = StructureContainer(cs)
for structure in structures:
    sc.add_structure(structure)

with open('log_files/3-structure-container', 'w') as f: # save the information in a file
    print('', sc, file=f)

# train the parameters, if necessary change the optimizer method
# using the 'fit_method' keyword in the Optimizer object
opt = Optimizer(sc.get_fit_data(), fit_method=optimization_method)
opt.train()

with open('log_files/4-optimizer', 'w') as f: # save the information in a file
    print('', opt, file=f)

# generate the force constant potential and save it in a .fcp file
fcp = ForceConstantPotential(cs, opt.parameters)
fcp.write('material.fcp')

with open('log_files/5-force-constant-potential', 'w') as f: # save the information in a file
    print('', fcp, file=f)
