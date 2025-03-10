# Pol Benítez Colominas, March 2025
# Universitat Politècnica de Catalunya

# Check if a phonon mode is polar or not

import yaml

from pymatgen.core import Structure
from pymatgen.analysis.bond_valence import BVAnalyzer

def get_oxidation_states(path):
    """
    It determines the oxidation states of all the atoms in a POSCAR file

    Inputs:
        path -> path to the POSCAR file
    """

    structure = Structure.from_file("POSCAR")
    
    bva = BVAnalyzer()
    structure_with_oxidation = bva.get_oxi_state_decorated_structure(structure)

    oxidations = []

    for site in structure_with_oxidation:
        for species in site.species.keys():
            charge = species.oxi_state
            oxidations.append(charge)

    return oxidations

def get_eigen_gamma(path):
    """
    Reads band.yaml file and return all the eigenvectors and eigenvalues of the phonon mode at gamma

    Inputs:
        path -> path to band.yaml file
    """

    with open(path, "r") as file:
        data = yaml.safe_load(file)

    number_phonons = data['natom'] * 3

    eigenvalues = []
    for mode in range(number_phonons):
        eigenvalues.append(data['phonon'][0]['band'][mode]['frequency'])

    eigenvectors = []
    for mode in range(number_phonons):
        eigenvector_atom = []
        for atom in range(number_phonons // 3):
            vec_x = data['phonon'][0]['band'][mode]['eigenvector'][atom][0][0]
            vec_y = data['phonon'][0]['band'][mode]['eigenvector'][atom][1][0]
            vec_z = data['phonon'][0]['band'][mode]['eigenvector'][atom][2][0]

            eigenvector_atom.append([vec_x, vec_y, vec_z])

        eigenvectors.append(eigenvector_atom)

    return eigenvalues, eigenvectors

def determine_polarity(oxidation, eigenvectors, threshold=1e-4):

    """
    It determines if a given phonon mode is polar or not (with a threshold accuracy)

    Inputs:
        oxidation -> list with the charge of all the atoms
        eigenvectors -> phonon mode vectors
        threshold -> accuracy of the polarity determination
   """

    component_x = 0
    component_y = 0
    component_z = 0

    for atom in range(len(oxidation)):
        component_x = component_x + (oxidation[atom] * eigenvectors[atom][0])
        component_y = component_y + (oxidation[atom] * eigenvectors[atom][1])
        component_z = component_z + (oxidation[atom] * eigenvectors[atom][2])

    if (abs(component_x) < threshold) and (abs(component_y) < threshold) and(abs(component_z) < threshold):
       return False
    else:
       return True

oxidations = get_oxidation_states('POSCAR')
eigenvalues, eigenvectors = get_eigen_gamma('band.yaml')

for mode in range(len(eigenvectors)):
    polar = determine_polarity(oxidations, eigenvectors[mode])
    print(f'The phonon mode {mode + 1} is polar: {polar}')
