# Pol Benítez Colominas, October 2024
# University of Cambridge and Universitat Politècnica de Catalunya

# Applies gaussian noise to distort a given structure and study the convergence with the functional

import os
import shutil

import numpy as np

from pymatgen.io.vasp import Poscar


#### Basic variables definition ####

width_disp = [0.50, 0.65, 0.80, 0.95, 1.10, 1.25, 1.40, 1.55, 1.70, 1.85] # standard deviation (width) of the gaussian noise
width_name = ['0_50', '0_65', '0_80', '0_95', '1_10', '1_25', '1_40', '1_55', '1_70', '1_85']
number_for_each = 4

path_to_save = 'generated-structures'    # path to save the resulting structures


#### Distorted structures generation ####
send_jobs = open(path_to_save + '/send_jobs.sh', "w")
send_jobs.write('#!/bin/bash \n')
send_jobs.write(' \n')

check_results = open(path_to_save + '/check_results.sh', "w")
check_results.write('#!/bin/bash \n')
check_results.write(' \n')

it_width = 0
poscar = Poscar.from_file('save/POSCAR')
original_poscar = poscar.structure
for disp in width_disp:
    for num in range(number_for_each):
        cp_poscar = original_poscar.copy()

        for site in cp_poscar:
            displacement_vector = np.random.uniform(0, disp, 3)

            site.coords = site.coords + displacement_vector

        os.mkdir(path_to_save + '/' + width_name[it_width] + '-' + str(num + 1).zfill(1))

        path_struc_PBEsol = path_to_save + '/' + width_name[it_width] + '-' + str(num + 1).zfill(1) + '/POSCAR'
        cp_poscar.to(filename=path_struc_PBEsol, fmt='Poscar')

        source = 'save/KPOINTS' 
        destination = path_to_save + '/' + width_name[it_width] + '-' + str(num + 1).zfill(1) + '/KPOINTS'
        shutil.copy(source, destination)

        source = 'save/INCAR' 
        destination = path_to_save + '/' + width_name[it_width] + '-' + str(num + 1).zfill(1) + '/INCAR'
        shutil.copy(source, destination)

        source = 'save/POTCAR' 
        destination = path_to_save + '/' + width_name[it_width] + '-' + str(num + 1).zfill(1) + '/POTCAR'
        shutil.copy(source, destination)

        source = 'save/run.sh' 
        destination = path_to_save + '/' + width_name[it_width] + '-' + str(num + 1).zfill(1) + '/run.sh'
        shutil.copy(source, destination)

        send_jobs.write('cd ' + width_name[it_width] + '-' + str(num + 1).zfill(1) + ' \n')
        send_jobs.write('sbatch run.sh \n')
        send_jobs.write('cd .. \n')
        send_jobs.write('\n')

        check_results.write('cd ' + width_name[it_width] + '-' + str(num + 1).zfill(1) + ' \n')
        check_results.write('echo ' + width_name[it_width] + '-' + str(num + 1).zfill(1) + '\n')
        check_results.write('tail -n 1 OSZICAR \n')
        check_results.write('cd .. \n')
        check_results.write('\n')

    it_width = it_width + 1

os.system("chmod +x " + path_to_save + "/send_jobs.sh")
os.system("chmod +x " + path_to_save + "/check_results.sh") 