# get materials centrosymmetric with large imaginary phonon

import os
import yaml

from pymatgen.io.vasp import Poscar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


def check_materials(structure, path):
    """
    Check if a material from Togo phonon data set is interesting or not
    Here we are interested in two conditions:
        - Centrosymmetry
        - Large imaginary phonon

    Inputs:
        structure: structure from Togo dataset to study
        path: path of the structure in Togo dataset
    """

    num_atoms = structure.num_sites
    # accept the material only if it has less than 10 atoms in the unit cell
    if num_atoms <= 10:

        # check if the space group contains inversion symmetry
        spacegroup_analyzer = SpacegroupAnalyzer(structure)
        centrosymmetric = spacegroup_analyzer.is_laue()

        # check if there is a large imaginary phonon
        with open(path + 'band.yaml', "r") as file:
            data = yaml.safe_load(file)

        symmetry_point = data['phonon'][0]['q-position']
        freqs = data['phonon'][0]['band'][0]['frequency']

        finish_condition = False
        it_q_point = 0
        while finish_condition == False:
            symmetry_point = data['phonon'][it_q_point]['q-position']

            if (symmetry_point[0] == 0) and (symmetry_point[1] == 0) and (symmetry_point[2] == 0):
                finish_condition = True
                freqs = data['phonon'][it_q_point]['band'][0]['frequency']
            else:
                it_q_point = it_q_point + 1
                freqs = 0

        if freqs <= -1.5:
            imaginary_phonon = True
        else:
            imaginary_phonon = False
    

        if (centrosymmetric == True) and (imaginary_phonon == True):
            valid_structure = True
        else:
            valid_structure = False

    else:
        valid_structure = False
        freqs = 0

    

    return valid_structure, freqs, num_atoms

directory_path = '../Loaded_PHONON'
folder_names = []

for folder_name in os.listdir(directory_path):
    folder_names.append(folder_name)

valid_structures = open('valid_structures1.txt', 'w')
valid_structures.write(f'mp-id     #atoms     frequency\n')

structures_names = []
for path_mat in folder_names:
    parts = path_mat.split('-')
    name_mp = '-'.join(parts[:3])
    structures_names.append(name_mp)


for it_struc in range(len(structures_names)):
    print(f'Structure number {it_struc + 1} of a total of {len(structures_names)} structures')
    path = directory_path + '/' + structures_names[it_struc] + '/'

    structure = Poscar.from_file(path + 'POSCAR').structure

    valid, freq, num = check_materials(structure, path)

    if valid == True:
        valid_structures.write(f'{structures_names[it_struc]}     {num}     {freq}\n')

valid_structures.close()