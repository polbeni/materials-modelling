# Pol Benítez Colominas, January 2024
# Universitat Politècnica de Catalunya

# Code to generate the force constants for different parameters (cutoffs and number of structures)

import os
import sys
import shutil

from ase.io import read, write
from hiphive import ClusterSpace, StructureContainer, ForceConstantPotential
from hiphive.utilities import prepare_structures
from trainstation import Optimizer

def get_directories_sorted(path, starts_with):
    """
    Creates an array with the names of the directories

    Inputs:
        path: path of the directories
        starts_with: common prefix in rattled structure directories

    Outputs:
        sorted_directories: array with the sorted directories
    """

    directories = [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d)) and d.startswith(starts_with)]
    
    def extract_number(directory):
        return int(directory.split('-')[-1]) 

    sorted_directories = sorted(directories, key=extract_number)
    
    return sorted_directories

def generate_force_constants(cutoffs, num_structures, optimization_method, data_path, results_path):
    """
    This function genereates the force constants for a number of structures and a given cutoffs

    Inputs: 
        cutoffs: cutoffs for pairs, triplets and quadruplets (in angstrom)
        num_structures: desired number of structures we want to take into account
        optimization_method: name of the desired optimization method employed
        data_path: distorted structures results path
        results_path: path where we want to save the results
    
    Outputs: 
        -
    """

    # read the structure input files
    prim = read('prim.extxyz') # primitive cell
    atoms_ideal = read('supercell_ideal.extxyz') # ideal supercell

    # create the cluster space for the defined cutoffs
    cs = ClusterSpace(prim, cutoffs)

    log_files_path = results_path + '/log_files'
    if os.path.exists(log_files_path):
        shutil.rmtree(log_files_path)
    os.mkdir(log_files_path)

    with open(log_files_path + '/1-cluster-space', 'w') as f: # save the information in a file
        print('', cs, file=f)

    original_stdout = sys.stdout
    with open(log_files_path + '/2-cluster-space-orbits', 'w') as file: # save the information in a file
        sys.stdout = file
        cs.print_orbits()
    sys.stdout = original_stdout

    # generate a list with the name of these directories
    prefix = 'disp'

    disp_direct = get_directories_sorted(data_path, prefix)

    # generate the rattled file, i.e., a file that contains energy, positions and
    # forces for each of the rattled structures
    rattled_file = open(results_path + '/rattled_structures.extxyz', "w")

    for disp in range(num_structures):
        path_structure = data_path + '/' + disp_direct[disp] + '/vasprun.xml'
        rattled = read(path_structure, format='vasp-xml')

        write(results_path + '/actual_rattled.extxyz', rattled)

        with open(results_path + '/actual_rattled.extxyz', 'r') as file:
            for line in file:
                rattled_file.write(line)

    os.remove(results_path + '/actual_rattled.extxyz')
    rattled_file.close()

    # read the rattled strucutres file
    rattled_structures = read(results_path + '/rattled_structures.extxyz', index=':')

    # generate the structure container
    structures = prepare_structures(rattled_structures, atoms_ideal)
    sc = StructureContainer(cs)
    for structure in structures:
        sc.add_structure(structure)

    with open(log_files_path + '/3-structure-container', 'w') as f: # save the information in a file
        print('', sc, file=f)

    # using the 'fit_method' keyword in the Optimizer object
    opt = Optimizer(sc.get_fit_data(), fit_method=optimization_method)
    opt.train()

    with open(log_files_path + '/4-optimizer', 'w') as f: # save the information in a file
        print('', opt, file=f)

    # generate the force constant potential and save it in a .fcp file
    fcp = ForceConstantPotential(cs, opt.parameters)
    fcp.write(results_path + '/material.fcp')

    with open(log_files_path + '/5-force-constant-potential', 'w') as f: # save the information in a file
        print('', fcp, file=f)

    return

# define all the cutoffs of interest and number of structures
diff_cutoffs = [[4, 4, 4], [4.5, 4, 4], [5, 4, 4]]
diff_num_structures = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500]

# define an array with all the directories names
dir_results = ['1_050', '1_100', '1_150', '1_200', '1_250', '1_300', '1_350', '1_400', '1_450', '1_500',
               '2_050', '2_100', '2_150', '2_200', '2_250', '2_300', '2_350', '2_400', '2_450', '2_500',
               '3_050', '3_100', '3_150', '3_200', '3_250', '3_300', '3_350', '3_400', '3_450', '3_500']

# define the loop to call generate_force_constants
num_iteration = 0
optm = 'least-squares'   # name of the desired optimization method employed
struc_path = 'path of the rattled structures results' # change this

for cuts in diff_cutoffs:
    for struc in diff_num_structures:
        print('########################################################################')
        print(f'Cutoffs: {cuts}')
        print(f'Number of structures: {struc}')

        if os.path.exists(dir_results[num_iteration]):
            shutil.rmtree(dir_results[num_iteration])
        os.mkdir(dir_results[num_iteration])

        generate_force_constants(cuts, struc, optm, struc_path, dir_results[num_iteration])

        num_iteration = num_iteration + 1

        print('########################################################################')
        print('')
        print('')
        print('')