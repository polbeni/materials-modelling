# Pol Benítez Colominas, March 2024
# Universitat Politècnica de Catalunya

# generates directories to study the energy convergence with different KPOINTS and ENCUT values
# it is necessary to provide the POSCAR, POTCAR and run.sh in a save directory

import os
import shutil


def generate_KPOINTS(path, kpoint):
    """
    Generates a KPOINT file centered in gamma with the specified number of kpoints

    Inputs:
        path -> path where we want to save the KPOINT file
        kpoint -> number of kpoints (a cubic kpoints grid will be generated)
    """

    KPOINTS = open(path + 'KPOINTS', 'w')
    KPOINTS.write('Automatic mesh \n')
    KPOINTS.write('0 \n')
    KPOINTS.write('G \n')
    KPOINTS.write(f' {int(kpoint)}   {int(kpoint)}   {int(kpoint)} \n')
    KPOINTS.write(' 0.  0.  0. \n')
    KPOINTS.close()

    return

def generate_INCAR(path, encut):
    """
    Generates a INCAR file with the specified ENCUT parameter

    Inputs:
        path -> path where we want to save the KPOINT file
        encut -> desired ENCUT parameter
    """

    INCAR = open(path + 'INCAR', 'w')
    INCAR.write('LCHARG   =  .FALSE. \n')
    INCAR.write('LWAVE    =  .FALSE. \n')
    INCAR.write('EDIFF    =  1E-7 \n')
    INCAR.write(f'ENCUT    =  {encut:.1f} \n')
    INCAR.write('ISMEAR   =  0 \n')
    INCAR.write('SIGMA    =  0.2 \n')
    INCAR.write('IALGO    =  48 \n')
    INCAR.write('NELM     =  200 \n')
    INCAR.write('GGA      =  PS \n')
    INCAR.close()

    return

possible_kpoints = [4, 5, 6, 7, 8, 9, 10] # desired kpoints grid
possible_encuts = [450, 500, 550, 600, 650, 700, 750] # desired encut values

path_save = 'save/'

send_jobs = open('send_jobs.sh', "w")
send_jobs.write('#!/bin/bash \n')
send_jobs.write(' \n')

check_results = open('check_results.sh', "w")
check_results.write('#!/bin/bash \n')
check_results.write(' \n')

for x in possible_kpoints:
    for y in possible_encuts:
        name_path = f'{x:01d}x{x:01d}x{x:01d}-{y:01d}'

        os.mkdir(name_path)

        generate_KPOINTS(name_path+'/', x)

        generate_INCAR(name_path+'/', y)

        source_file = path_save + '/POSCAR'
        destination_file = name_path + '/POSCAR'
        shutil.copy(source_file, destination_file)

        source_file = path_save + '/POTCAR'
        destination_file = name_path + '/POTCAR'
        shutil.copy(source_file, destination_file)

        source_file = path_save + '/run.sh'
        destination_file = name_path + '/run.sh'
        shutil.copy(source_file, destination_file)

        send_jobs.write('cd ' + name_path + ' \n')
        send_jobs.write('sbatch run.sh \n')
        send_jobs.write('cd .. \n')
        send_jobs.write('\n')

        check_results.write('cd ' + name_path + ' \n')
        check_results.write('tail -n 1 OUTCAR \n')
        check_results.write('cd .. \n')
        check_results.write('\n')

os.system("chmod +x send_jobs.sh")
os.system("chmod +x check_results.sh")