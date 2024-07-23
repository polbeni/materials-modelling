import os
import shutil

from pymatgen.io.vasp import Poscar

def generate_POSCAR(path, stoi, path_to_save):
    """
    Generates a POSCAR with the desired VCA

    Inputs:
        path: path to the POSCAR file
        stoi: desired stoichiometry for the VCA
        path_to_save: path to save the new POSCAR generated
    """

    poscar_file = open(path, 'r')
    new_poscar_file = open(path_to_save, 'w')

    poscar_file.readline()
    new_poscar_file.write('VCA POSCAR file\n')

    for x in range(4):
        line = poscar_file.readline()
        new_poscar_file.write(line)

    atoms_line = 'Ag  '
    atoms_number = '3   '

    if stoi[0] != 0:
        atoms_line = atoms_line + 'S  '
        atoms_number = atoms_number + '1  '

    if stoi[1] != 0:
        atoms_line = atoms_line + 'Se  '
        atoms_number = atoms_number + '1  '

    if stoi[2] != 0:
        atoms_line = atoms_line + 'Br  '
        atoms_number = atoms_number + '1  '

    if stoi[3] != 0:
        atoms_line = atoms_line + 'I  '
        atoms_number = atoms_number + '1  '

    new_poscar_file.write(atoms_line + '\n')
    new_poscar_file.write(atoms_number + '\n')

    for x in range(2):
        poscar_file.readline()

    for x in range(4):
        line = poscar_file.readline()
        new_poscar_file.write(line)

    if stoi[0] != 0:
        new_poscar_file.write('  0.5000000000000000  0.5000000000000000  0.5000000000000000\n')

    if stoi[1] != 0:
        new_poscar_file.write('  0.5000000000000000  0.5000000000000000  0.5000000000000000\n')

    if stoi[2] != 0:
        new_poscar_file.write(' -0.0000000000000000 -0.0000000000000000  0.0000000000000000\n')

    if stoi[3] != 0:
        new_poscar_file.write(' -0.0000000000000000 -0.0000000000000000  0.0000000000000000\n')

    return

def generate_INCAR(path, stoi, path_to_save):
    """
    Generates a INCAR with the desired VCA flag

    Inputs:
        path: path to the INCAR file
        stoi: desired stoichiometry for the VCA
        path_to_save: path to save the new INCAR generated
    """

    new_incar_file = open(path_to_save, 'w')

    with open(path, 'r') as file:
        for line in file:
            new_incar_file.write(line)

    vca_line = 'VCA      =  '
    vca_line = vca_line + '1.0 '

    if stoi[0] != 0:
        vca_line = vca_line + str(format(stoi[0], '.2f')) + ' '

    if stoi[1] != 0:
        vca_line = vca_line + str(format(stoi[1], '.2f')) + ' '

    if stoi[2] != 0:
        vca_line = vca_line + str(format(stoi[2], '.2f')) + ' '

    if stoi[3] != 0:
        vca_line = vca_line + str(format(stoi[3], '.2f')) + ' '

    new_incar_file.write(vca_line)

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
    gga_path = '/home/pol/work/materials-modelling/pseudopotentials-files/VASP/GGA/'

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

S_coef = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
Br_coef = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]

if os.path.exists('VCA_structures'):
    shutil.rmtree('VCA_structures')
os.mkdir('VCA_structures')

send_jobs = open('VCA_structures/send_jobs-relax.sh', "w")
send_jobs.write('#!/bin/bash \n')
send_jobs.write(' \n')

check_results = open('VCA_structures/check_results-relax.sh', "w")
check_results.write('#!/bin/bash \n')
check_results.write(' \n')

iteration = 1
for S_in in S_coef:
    for Br_in in Br_coef:
        stoi_array = [S_in, 1-S_in, Br_in, 1-Br_in]

        os.mkdir('VCA_structures/vca-' + str(iteration).zfill(3))
        path_sim = 'VCA_structures/vca-' + str(iteration).zfill(3) + '/'

        ### Relaxation 
        os.mkdir(path_sim + 'relaxation')
        path_relax = path_sim + 'relaxation/'

        ### Band gap 
        os.mkdir(path_sim + 'bandgap')
        path_bg = path_sim + 'bandgap/'

        # POSCAR
        generate_POSCAR('POSCAR', stoi_array, path_relax + '/POSCAR')
        
        # INCAR
        generate_INCAR('save-relaxation/INCAR', stoi_array, path_relax + '/INCAR')
        generate_INCAR('save-bandgap/INCAR', stoi_array, path_bg + '/INCAR')

        # KPOINTS
        source = 'save-relaxation/KPOINTS'
        destination = path_relax + '/KPOINTS'
        shutil.copy(source, destination)
        source = 'save-bandgap/KPOINTS'
        destination = path_bg + '/KPOINTS'
        shutil.copy(source, destination)

        # run.sh
        source = 'save-relaxation/run.sh'
        destination = path_relax + '/run.sh'
        shutil.copy(source, destination)
        source = 'save-bandgap/run.sh'
        destination = path_bg + '/run.sh'
        shutil.copy(source, destination)

        # POTCAR
        atoms_list = get_atoms_array(path_relax + '/POSCAR')
        create_potcar(atoms_list, path_relax + '/POTCAR')
        create_potcar(atoms_list, path_bg + '/POTCAR')

        # scripts
        send_jobs.write('cd ' + 'vca-' + str(iteration).zfill(3) + '/relaxation \n')
        send_jobs.write('sbatch run.sh \n')
        send_jobs.write('cd ../.. \n')
        send_jobs.write('\n')

        check_results.write('cd ' + 'vca-' + str(iteration).zfill(3) + '/relaxation \n')
        check_results.write('echo ' 'vca-' + str(iteration).zfill(3) + '\n')
        check_results.write('tail -n 1 OSZICAR \n')
        check_results.write('cd ../.. \n')
        check_results.write('\n')

        iteration = iteration + 1

os.system("chmod +x VCA_structures/send_jobs-relax.sh")
os.system("chmod +x VCA_structures/check_results-relax.sh")