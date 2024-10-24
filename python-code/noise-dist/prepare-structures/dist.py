# Pol Benítez Colominas, October 2024
# University of Cambridge and Universitat Politècnica de Catalunya

# Applies gaussian noise to distort a given structure and study the convergence with the functional

import os
import shutil

import numpy as np

from pymatgen.io.vasp import Poscar


#### Basic variables definition ####

material = ['Ag3SBr', 'Ag3SI']
supercell = ['1x1x1', '2x2x2']
width_disp = [0.5, 0.9, 1.3, 1.7] # standard deviation (width) of the gaussian noise
number_for_each = 5

path_to_save = 'generated-structures'    # path to save the resulting structures


#### Distorted structures generation ####
send_jobs = open(path_to_save + '/send_jobs.sh', "w")
send_jobs.write('#!/bin/bash \n')
send_jobs.write(' \n')

check_results = open(path_to_save + '/check_results.sh', "w")
check_results.write('#!/bin/bash \n')
check_results.write(' \n')

num_struc = 1
for mat in material:
    for cell in supercell:
        poscar_path = 'original-data/' + mat + '/POSCAR-' + cell
        poscar = Poscar.from_file(poscar_path)
        original_poscar = poscar.structure
        for disp in width_disp:
            for num in range(number_for_each):
                cp_poscar = original_poscar.copy()

                for site in cp_poscar:
                    displacement_vector = np.random.uniform(0, disp, 3)

                    site.coords = site.coords + displacement_vector

                os.mkdir(path_to_save + '/struc-' + str(num_struc).zfill(3))
                os.mkdir(path_to_save + '/struc-' + str(num_struc).zfill(3) + '/PBEsol')
                os.mkdir(path_to_save + '/struc-' + str(num_struc).zfill(3) + '/HSEsol')

                path_struc_PBEsol = path_to_save + '/struc-' + str(num_struc).zfill(3) + '/PBEsol/POSCAR'
                cp_poscar.to(filename=path_struc_PBEsol, fmt='Poscar')
                path_struc_HSEsol = path_to_save + '/struc-' + str(num_struc).zfill(3) + '/HSEsol/POSCAR'
                cp_poscar.to(filename=path_struc_HSEsol, fmt='Poscar')

                source = 'save/KPOINTs/s' + cell
                destination = path_to_save + '/struc-' + str(num_struc).zfill(3) + '/PBEsol/KPOINTS'
                shutil.copy(source, destination)
                destination = path_to_save + '/struc-' + str(num_struc).zfill(3) + '/HSEsol/KPOINTS'
                shutil.copy(source, destination)

                source = 'save/POTCARs/' + mat
                destination = path_to_save + '/struc-' + str(num_struc).zfill(3) + '/PBEsol/POTCAR'
                shutil.copy(source, destination)
                destination = path_to_save + '/struc-' + str(num_struc).zfill(3) + '/HSEsol/POTCAR'
                shutil.copy(source, destination)

                source = 'save/PBEsol/INCAR'
                destination = path_to_save + '/struc-' + str(num_struc).zfill(3) + '/PBEsol/INCAR'
                shutil.copy(source, destination)
                source = 'save/PBEsol/run.sh'
                destination = path_to_save + '/struc-' + str(num_struc).zfill(3) + '/PBEsol/run.sh'
                shutil.copy(source, destination)

                source = 'save/HSEsol/INCAR'
                destination = path_to_save + '/struc-' + str(num_struc).zfill(3) + '/HSEsol/INCAR'
                shutil.copy(source, destination)
                source = 'save/HSEsol/run.sh'
                destination = path_to_save + '/struc-' + str(num_struc).zfill(3) + '/HSEsol/run.sh'
                shutil.copy(source, destination)

                num_struc = num_struc + 1

                send_jobs.write('cd ' + '/struc-' + str(num_struc).zfill(3) + '/PBEsol' + ' \n')
                send_jobs.write('sbatch run.sh \n')
                send_jobs.write('cd ../HSEsol \n')
                send_jobs.write('sbatch run.sh \n')
                send_jobs.write('cd ../.. \n')
                send_jobs.write('\n')

                check_results.write('cd ' + '/struc-' + str(num_struc).zfill(3) + '/PBEsol' + ' \n')
                check_results.write('echo ' + 'struc-' + str(num_struc).zfill(3) + ' PBEsol' + '\n')
                check_results.write('tail -n 1 OSZICAR \n')
                check_results.write('cd ../HSEsol \n')
                check_results.write('echo ' + 'struc-' + str(num_struc).zfill(3) + ' HSEsol' + '\n')
                check_results.write('tail -n 1 OSZICAR \n')
                check_results.write('cd ../.. \n')
                check_results.write('\n')

os.system("chmod +x " + path_to_save + "/send_jobs.sh")
os.system("chmod +x " + path_to_save + "/check_results.sh") 