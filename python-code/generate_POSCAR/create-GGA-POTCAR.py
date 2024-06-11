# Pol Benítez Colominas, June 2024
# Universitat Politècnica de Catalunya

# This code generates a POTCAR file for a given POSCAR file

import os
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

atoms_list = get_atoms_array('materials/mp-626151')
create_potcar(atoms_list, 'POTCAR')