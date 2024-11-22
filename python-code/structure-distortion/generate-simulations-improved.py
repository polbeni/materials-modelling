# Pol Benítez Colominas, September 2024
# Universitat Politècnica de Catalunya and University of Cambridge

# This code generates the distorted POSCAR files for the desired eigenvectors at gamma-point in
# order to study LO-TO splitting effect in e-p band gap

import os
import shutil
import subprocess

def execute_claudi(path, amplitude):
    """
    Executes Claudi program to distort the structure

    Inputs:
        path: path where the calculation is performed
    """
    #exe_path = './displ.exe'

    input1 = '3' # number of different species in the unit cell
    input2 = '5' # total number of atoms in the unit cell
    input3 = amplitude # amplitude displacment in angstroms
    input4 = '0' # to include normalization of eigenvectors or not (0: No, 1: Yes)
    inputs = f'{input1}\n{input2}\n{input3}\n{input4}\n'

    result = subprocess.run([path], input=inputs, text=True, capture_output=True)

    return

def generate_VECTORS(path_OUTCAR, num_phonon, output_path):
    """
    Generates a VECTORS file (with Claudi's script format) reading the eigenvalues of the phonon modes at gamma-point 
    from the OUTCAR file for the desired phonon mode

    Inputs:
        path_OUTCAR -> path to the OUTCAR file
        num_phonon -> number of the phonon mode
        output_path -> path where to save VECTORS file
    """

    with open(path_OUTCAR, 'r') as file:
        searc_pattern = 'Eigenvectors and eigenvalues of the dynamical matrix'

        num_line = 1
        for line in file:
            
            if searc_pattern in line:
                val_line = num_line
            else:
                num_line = num_line + 1

    OUTCAR = open(path_OUTCAR, 'r')
    for _ in range(val_line):
        OUTCAR.readline()
    
    for _ in range(3):
        OUTCAR.readline()

    for _ in range(8 * (num_phonon - 1)):
        OUTCAR.readline()

    for _ in range(2):
        OUTCAR.readline()

    VECTORS = open(output_path, 'w')
    
    for x in range(5): # iterate over all the atoms in the unit cell
        line = OUTCAR.readline()
        VECTORS.write(line)

    VECTORS.close()
    OUTCAR.close()
    
    return 

            

list_amplitudes = ['0.2', '0.4', '0.6']

phonons_calc = ['phonon-01', 'phonon-02', 'phonon-03', 'phonon-04', 'phonon-05', 
                'phonon-06', 'phonon-07', 'phonon-08', 'phonon-09', 'phonon-10', 
                'phonon-11', 'phonon-12', 'phonon-13', 'phonon-14', 'phonon-15']

ampl = ['2', '4', '6']


send_jobs = open('send_jobs.sh', "w")
send_jobs.write('#!/bin/bash \n')
send_jobs.write(' \n')

check_results = open('check_results.sh', "w")
check_results.write('#!/bin/bash \n')
check_results.write(' \n')

for phonon in range(15):
    for amp_value in range(3):

        generate_VECTORS('OUTCAR', phonon + 1, 'VECTORS')
        execute_claudi('./displ.exe', list_amplitudes[amp_value])

        dir_name = phonons_calc[phonon] + '-' + ampl[amp_value]
        print(dir_name)

        os.mkdir(dir_name)

        # POSCAR
        source = 'POSCARnew'
        destination = dir_name + '/POSCAR'
        shutil.copy(source, destination)

        # KPOINTS
        source = 'save/KPOINTS'
        destination = dir_name + '/KPOINTS'
        shutil.copy(source, destination)

        # INCAR
        source = 'save/INCAR'
        destination = dir_name + '/INCAR'
        shutil.copy(source, destination)

        # POTCAR
        source = 'save/POTCAR'
        destination = dir_name + '/POTCAR'
        shutil.copy(source, destination)

        # run.sh
        source = 'save/run.sh'
        destination = dir_name + '/run.sh'
        shutil.copy(source, destination)

        send_jobs.write('cd ' + dir_name + '\n')
        send_jobs.write('sbatch run.sh \n')
        send_jobs.write('cd .. \n')
        send_jobs.write('\n')

        check_results.write('cd ' + dir_name + '\n')
        check_results.write('echo ' + dir_name +  ' \n')
        check_results.write('tail -n1 OSZICAR \n')
        check_results.write('cd .. \n')
        check_results.write('\n')

send_jobs.close()
os.system("chmod +x " + "send_jobs.sh")

check_results.close()
os.system("chmod +x " + "check_results.sh")