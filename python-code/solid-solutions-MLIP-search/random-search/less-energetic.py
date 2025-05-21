# Pol Benítez Colominas, May 2025
# Universitat Politècnica de Catalunya

# Computes random solutions permuations and relaxes the structures
# and compute their energies with MLIPs (MACE-based)

import os
import shutil
import itertools
from copy import deepcopy
import random

import warnings

from pymatgen.core import Structure
from pymatgen.io.vasp import Poscar
from pymatgen.io.ase import AseAtomsAdaptor

from ase.optimize import BFGS, FIRE
from ase.io.trajectory           import Trajectory
from ase.io.vasp                 import write_vasp
from ase.constraints import ExpCellFilter, FixAtoms

from mace.calculators import mace_mp
import torch
torch.serialization.add_safe_globals([slice])

warnings.simplefilter('ignore')


######################## Variables #########################
atoms_solid_solution1 = 'Br'
atoms_solid_solution2 = 'I'
ss_coef_Br = 0.2
num_random_struc = 10
############################################################


######################## MACE stuff ########################
if torch.cuda.is_available():
    device = 'cuda'
    print("GPU is available. Using GPU.")
else:
    device = 'cpu'
    print("GPU not available. Using CPU.")

model = 'large'    # name of the model
fmax = 0.05        # max force allowed
smax = 200         # max number of steps for the relaxation
############################################################


######################### Functions ########################
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
############################################################


########################### Main ###########################
# Open the POSCAR
structure = Poscar.from_file('POSCAR').structure

list_without, list_with = get_atoms_list(structure, atoms_solid_solution1)

# Generate a path to save the results
if os.path.exists('results'):
    shutil.rmtree('results')
os.mkdir('results')

# Loop random search
list_struc_values = []

for it_struc in range(num_random_struc):
    # Generate the proportion of solid solution
    list_with_solid_solution = generate_list_solid_solution(list_with, atoms_solid_solution2, ss_coef_Br)

    # Shuffle the solid solution atoms in order to have a random solid solution supercell
    random.shuffle(list_with_solid_solution)

    # Generate the random solid solution
    solid_solution = deepcopy(structure)

    solid_solution_atoms = list_without + list_with_solid_solution

    new_species, new_coords = sort_atoms(solid_solution_atoms, solid_solution.frac_coords, atoms_solid_solution1, atoms_solid_solution2)

    solid_solution_structure = Structure(
        lattice=solid_solution.lattice,
        species=new_species,
        coords=new_coords,
        coords_are_cartesian=False
    )

    # Create ase object and perform md with MACE
    ase_adaptor = AseAtomsAdaptor()
    atoms = ase_adaptor.get_atoms(solid_solution_structure)

    atoms.calc = mace_mp(model=model, device=device)

    # Fix the symmetry (just change the lattice parameters)
    constraint = FixAtoms(indices=[atom.index for atom in atoms]) # fix the atoms positions in the relaxation
    atoms.set_constraint(constraint)

    mask = [1, 1, 1, 0, 0, 0]
    atoms_filter = ExpCellFilter(atoms, mask=mask) # allow lattice parameters to change

    # Perform the relaxation
    dyn = BFGS(atoms_filter) # relax the structure
    dyn.run(fmax=fmax, steps=smax)

    # Save the relaxed solid solution and fix the POSCAR format
    write_vasp('POSCAR-after-relax', atoms=atoms, direct=True)

    old_POSCAR = open('POSCAR-after-relax', 'r') # fix the POSCAR structure format
    new_POSCAR = open('results/' + 'POSCAR-' + str(it_struc + 1).zfill(3), 'w')

    for x in range(7):
        line = old_POSCAR.readline()
        new_POSCAR.write(line)

        if x == 6: # read the number of atoms
            num_atoms = 0
            for atom in range(len(line.split())):
                num_atoms = num_atoms + int(line.split()[atom])

    old_POSCAR.readline()

    line = old_POSCAR.readline()
    new_POSCAR.write(line)

    for x in range(num_atoms):
        line = old_POSCAR.readline()
        new_POSCAR.write(f'   {float(line.split()[0])}   {float(line.split()[1])}   {float(line.split()[2])}\n')

    old_POSCAR.close()
    new_POSCAR.close()

    # Determine the energy of the structure and save the value
    energy = atoms.get_potential_energy()

    list_struc_values.append([it_struc + 1, energy / len(solid_solution_atoms)])

sorted_list = sorted(list_struc_values, key=lambda x: x[1], reverse=False)

with open('results/' + 'results_file.txt', 'w') as file:
    file.write('Structure       E(eV/atom)\n')
    
    for ss in range(len(sorted_list)):
        file.write(f'{sorted_list[ss][0]}       {sorted_list[ss][1]}\n')
############################################################