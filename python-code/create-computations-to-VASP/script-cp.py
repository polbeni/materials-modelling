import os
import shutil
from pymatgen.io.vasp import Poscar

def get_atoms_array(poscar_file):
    """
    Returns an array with the unique atoms in POSCAR

    Inputs:
        poscar_file: path to the POSCAR file (and name of the file)
    """

    poscar = Poscar.from_file(poscar_file)

    element_symbols = poscar.site_symbols

    unique_elements = []
    for element in element_symbols:
        if element not in unique_elements:
            unique_elements.append(element)

    return unique_elements


def create_potcar(atoms, path_potcar):
    """
    Creates a POTCAR file (for now with GGA functionals) for the desired atoms

    Inputs:
        atoms: array with atoms
        path_potcar: path to save the POTCAR file (and name of the file)
    """
    gga_path = '/home/pol/Documents/work/materials-modelling/pseudopotentials-files/VASP/GGA/'

    potcar_generation_command = 'cat '

    for atom in atoms:
        if atom in ['Ba', 'Ca', 'Cs', 'Fr', 'Nb', 'Ra', 'Rb', 'Sr', 'Y']: # check exceptions
            atom_path = gga_path + atom + '_sv' + '/POTCAR'
        elif atom == 'K':
            atom_path = gga_path + atom + '_pv' + '/POTCAR'
        else:
            atom_path = gga_path + atom + '/POTCAR'

        potcar_generation_command = potcar_generation_command + atom_path + ' '

    potcar_generation_command = potcar_generation_command + '> ' + path_potcar

    os.system(potcar_generation_command)

    return

materials = ['Ag3SCl', 'Ag3SBr', 'Ag3SI', 'Ag3SeCl', 'Ag3SeBr', 'Ag3SeI']

for material in materials:
    path_struc = 'structures/' + material + '/'
    path_results = 'results/' + material + '/'

    send_jobs = open(path_results + 'send_jobs.sh', "w")
    send_jobs.write('#!/bin/bash \n')
    send_jobs.write(' \n')

    check_results = open(path_results + 'check_results.sh', "w")
    check_results.write('#!/bin/bash \n')
    check_results.write(' \n')
    
    for x in range(10):
        os.makedirs(path_results + 'simulation-' + str(x + 1).zfill(3))

        # POSCAR
        source = path_struc + 'POSCAR-' + str(x + 1).zfill(3)
        destination = path_results + 'simulation-' + str(x + 1).zfill(3) + '/POSCAR'
        shutil.copy(source, destination)

        # POTCAR
        atoms_list = get_atoms_array(path_results + 'simulation-' + str(x + 1).zfill(3) + '/POSCAR')
        create_potcar(atoms_list, path_results + 'simulation-' + str(x + 1).zfill(3) + '/POTCAR')

        # KPOINTS
        source = 'save/KPOINTS'
        destination = path_results + 'simulation-' + str(x + 1).zfill(3) + '/KPOINTS'
        shutil.copy(source, destination)

        # run.sh and INCAR
        if ((x + 2)%2) == 0:
            source = 'save/INCAR-bg'
            destination = path_results + 'simulation-' + str(x + 1).zfill(3) + '/INCAR'
            shutil.copy(source, destination)

            source = 'save/run-bg.sh'
            destination = path_results + 'simulation-' + str(x + 1).zfill(3) + '/run.sh'
            shutil.copy(source, destination)
        else:
            source = 'save/INCAR-opt'
            destination = path_results + 'simulation-' + str(x + 1).zfill(3) + '/INCAR'
            shutil.copy(source, destination)

            source = 'save/run-opt.sh'
            destination = path_results + 'simulation-' + str(x + 1).zfill(3) + '/run.sh'
            shutil.copy(source, destination)

        # scripts
        send_jobs.write('cd ' + 'simulation-' + str(x + 1).zfill(3) + ' \n')
        send_jobs.write('sbatch run.sh \n')
        send_jobs.write('cd .. \n')
        send_jobs.write('\n')

        check_results.write('cd ' + 'simulation-' + str(x + 1).zfill(3) + ' \n')
        check_results.write('echo ' + 'simulation-' + str(x + 1).zfill(3) + '\n')
        check_results.write('tail -n 1 OSZICAR \n')
        check_results.write('cd .. \n')
        check_results.write('\n')

    os.system("chmod +x " + path_results + "send_jobs.sh")
    os.system("chmod +x " + path_results + "check_results.sh")   