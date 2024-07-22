# Pol Benítez Colominas, July 2024
# Universitat Politècnica de Catalunya

# Generates all the possible permutations for a solid solution
# it does not take into account equivalent structures

import os
import shutil
import itertools
from copy import deepcopy

from pymatgen.io.vasp import Poscar
from pymatgen.core import Structure

def find_all_permutations(elements):
    """
    Find all the possible permutations for a list of elements and returns all the possible lists

    Inputs:
        elements: array with the elements
    """

    unique_permutations = set(itertools.permutations(elements))

    unique_permutations_list = [list(p) for p in unique_permutations]

    return unique_permutations_list

def get_atoms_list(stucture, atom_type):
    """
    Returns two list of atoms for a given structure, one with all the atoms != to atom, and another all
    the atoms == to atom

    Inputs:
        structure: unit cell structure in pymatgen format
        atom_type: the type of atom    
    """

    atom_list = [str(site.specie) for site in structure.sites]

    without_atom = []
    with_atom = []

    for atom in atom_list:
        if atom_type == atom:
            with_atom.append(atom)
        else:
            without_atom.append(atom)

    return without_atom, with_atom

def generate_list_solid_solution(list_atoms, atom_type, stoi_coef):
    """
    Generates a list with the atom in the solid solution with the desired stoichiometry coefficient

    Inputs:
        list_atoms: list of atoms to modify
        atom_type: the other type of atom of the solid solution
        stoi_coef: coefficient of the stoichiometry of the solid solution
    """

    new_list_atoms = []

    for x in range(len(list_atoms)):
        if (x + 1)/len(list_atoms) > stoi_coef:
            new_list_atoms.append(atom_type)
        else:
            new_list_atoms.append(list_atoms[x])

    return new_list_atoms

def sort_atoms(species, coords, atom1, atom2):
    """
    It recives a list of species and coords and returns both lists ordered first by atom1 and then atom2
    """

    new_species = []
    new_coords = []

    for x in range(len(species)):
        if (species[x] != atom1) and (species[x] != atom2):
            new_species.append(species[x])
            new_coords.append(coords[x])

    for x in range(len(species)):
        if species[x] == atom1:
            new_species.append(species[x])
            new_coords.append(coords[x])

    for x in range(len(species)):
        if species[x] == atom2:
            new_species.append(species[x])
            new_coords.append(coords[x])

    return new_species, new_coords

    
######### PARAMETERS #########
size_supercell = 2
atoms_solid_solution1 = 'Br'
atoms_solid_solution2 = 'I'
solid_solution_coef = 0.5
##############################


structure = Poscar.from_file('POSCAR').structure

scaling_matrix = [[size_supercell, 0, 0], [0, size_supercell, 0], [0, 0, size_supercell]]
supercell = structure.make_supercell(scaling_matrix)

list_without, list_with = get_atoms_list(supercell, atoms_solid_solution1)

list_with_solid_solution = generate_list_solid_solution(list_with, atoms_solid_solution2, solid_solution_coef)

unique_permutations_list = find_all_permutations(list_with_solid_solution)

print(f"Number of unique permutations: {len(unique_permutations_list)}")

if os.path.exists('solid_solution_structures'):
    shutil.rmtree('solid_solution_structures')
os.mkdir('solid_solution_structures')

iteration = 1
for permutation in unique_permutations_list:
    solid_solutions_atoms = list_without + permutation

    solid_solution = deepcopy(supercell)

    new_species, new_coords = sort_atoms(solid_solutions_atoms, solid_solution.frac_coords, atoms_solid_solution1, atoms_solid_solution2)
    
    solid_solution_structure = Structure(
        lattice=solid_solution.lattice,
        species=new_species,
        coords=new_coords,
        coords_are_cartesian=False
    )
    
    # Write the modified structure to a new POSCAR file
    solid_solution_structure = Poscar(solid_solution_structure)
    solid_solution_structure.write_file('solid_solution_structures/POSCAR-' + str(iteration).zfill(3))

    iteration = iteration + 1