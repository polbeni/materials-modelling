# Pol Benítez Colominas, November 2024
# University of Cambridge and Universitat Politècnica de Catalunya

# Reads a phonopy.yaml file and returns the original structure in POSCAR file format

import yaml
from collections import Counter

def get_POSCAR(path_phonopy, path_final):
    """
    It generates a POSCAR of the original structure from the phonopy.yaml file

    Inputs:
        path_phonopy: path the the phonopy.yaml file
        path_final: path and name of the resulting POSCAR file
    """

    # load the phonopy.yaml file
    with open(path_phonopy, 'r') as f:
        data = yaml.safe_load(f)

    # first two lines
    POSCAR = open(path_final, 'w')
    POSCAR.write('POSCAR from phonopy.yaml file \n')
    POSCAR.write('1.0 \n')

    # lattice parameters
    lattice = data['primitive_cell']['lattice']
    POSCAR.write(f'   {lattice[0][0]:.16f}   {lattice[0][1]:.16f}   {lattice[0][2]:.16f} \n')
    POSCAR.write(f'   {lattice[1][0]:.16f}   {lattice[1][1]:.16f}   {lattice[1][2]:.16f} \n')
    POSCAR.write(f'   {lattice[2][0]:.16f}   {lattice[2][1]:.16f}   {lattice[2][2]:.16f} \n')

    # atoms type and number
    atom_info = data['primitive_cell']['points']
    atom_symbols = [atom['symbol'] for atom in atom_info]
    atom_counts = Counter(atom_symbols)

    list_atoms = []
    list_number = []
    for atom_type, count in atom_counts.items():
        list_atoms.append(atom_type)
        list_number.append(count)

    line_atom = ''
    for atom in list_atoms:
        line_atom = line_atom + atom + ' '

    line_number = ''
    for num in list_number:
        line_number = line_number + str(num) + ' '
    
    POSCAR.write(f'{line_atom} \n')
    POSCAR.write(f'{line_number} \n')

    # Direct line
    POSCAR.write('Direct \n')

    # Atoms positions
    total_num_atoms = 0
    for num in list_number:
        total_num_atoms = total_num_atoms + num

    for atom in range(total_num_atoms):
        positions = data['primitive_cell']['points'][atom]['coordinates']
        POSCAR.write(f'   {positions[0]:.16f}   {positions[1]:.16f}   {positions[2]:.16f} \n')

    return


get_POSCAR('phonopy.yaml', 'POSCAR')