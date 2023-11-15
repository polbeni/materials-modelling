# Pol Benítez Colominas, November 2023
# Universitat Politècnica de Catalunya

# Code to generate rattled POSCAR files to send to VASP with hiPhive
# when cell_size=1, there is no bug in ion sorting

import os

from ase.io import write
from ase.io.vasp import read_vasp
from hiphive.structure_generation import generate_mc_rattled_structures

# parameters
n_structures = 30 # number of structures
cell_size = 1 # size of the supercell
rattle_std = 0.03 # standard deviation of the distribution of displacements
minimum_distance = 2.3 # minimum separation between two atoms in the rattled structures

# read POSCAR file
prim = read_vasp(file='POSCAR')
atoms_ideal = prim.repeat(cell_size)

# generate the structures
structures = generate_mc_rattled_structures(atoms_ideal, n_structures, rattle_std, minimum_distance)

# save the structures for later use in force determination
write('prim.extxyz', prim) # primitive cell
write('supercell_ideal.extxyz', atoms_ideal) # ideal supercell
# write('supercells_rattled.extxyz', structures)

"""
# fix the ion sort problem in atoms_ideal file
supercell_corrected = open('supercell_ideal.extxyz', 'w')
supercell_prev = open('supercell_ideal_prev.extxyz', 'r')

for x in range(2):
    line_prev = supercell_prev.readline()
    supercell_corrected.write(line_prev)

    if x == 0:
        total_atoms = int(line_prev)

supercell_prev.close()

poscar = open('POSCAR', 'r')
for x in range(5):
    poscar.readline()
atoms_line = poscar.readline()
atoms = list(atoms_line.split())
stoichiometry_line = poscar.readline()
stoichiometry = list(stoichiometry_line.split())
poscar.close()

total_num_atoms_stoi = 0
for stoi in stoichiometry:
    total_num_atoms_stoi = total_num_atoms_stoi + int(stoi)

num_iteration = 1
for actual_stoi in stoichiometry:
    supercell_prev = open('supercell_ideal_prev.extxyz', 'r')
    for x in range(2):
        supercell_prev.readline()

    if num_iteration == 1:
        previous_stoi = 0
        next_stoi = total_num_atoms_stoi - int(actual_stoi)
    elif num_iteration == len(stoichiometry):
        previous_stoi = total_num_atoms_stoi - int(actual_stoi)
        next_stoi = 0
    else:
        previous_stoi = 0
        for y in range(num_iteration - 1):
            previous_stoi = previous_stoi + int(stoichiometry[y])
        next_stoi = total_num_atoms_stoi - int(actual_stoi) - previous_stoi

    atoms_considered = 0
    while atoms_considered <= total_atoms:
        if previous_stoi != 0:
            for y in range(previous_stoi):
                supercell_prev.readline()
                atoms_considered = atoms_considered + 1

        for y in range(int(actual_stoi)):
            atom_position_line = supercell_prev.readline()
            supercell_corrected.write(atom_position_line)
            atoms_considered = atoms_considered + 1

        if next_stoi != 0:
            for y in range(next_stoi):
                supercell_prev.readline()
                atoms_considered = atoms_considered + 1

    num_iteration = num_iteration + 1


supercell_prev.close()
supercell_corrected.close()

os.remove('supercell_ideal_prev.extxyz')
"""

# save the structures in POSCAR file (some corrections have to be done)
for x in range(n_structures):
    structures[x].write(f'POSCAR-{x+1:04}')

"""
# fix the POSCAR files
for x in range(n_structures):
    # open files
    pre_file = open(f'POSCAR-{x+1:04}_pre', "r")
    final_file = open(f'POSCAR-{x+1:04}', "w")

    # header
    pre_file.readline()
    final_file.write(f'Generated structure number {x+1:04} with Hiphive\n')

    # scaling factor and lattice
    for x_line in range(4):
        line_pre = pre_file.readline()
        final_file.write(line_pre)

    # atoms
    line_pre = pre_file.readline()
    atoms_line = line_pre.split()
    print(atoms_line)
    is_it_condition = False
    diff_atoms = []
    counter_element = 0
    while is_it_condition == False:
        print(counter_element)
        actual_element = atoms_line[counter_element]
        if actual_element in diff_atoms:
            is_it_condition = True
            break
        diff_atoms.append(actual_element)
        counter_element = counter_element + 1
    num_diff_atoms = len(diff_atoms)
    
    new_line = ''
    for atom in diff_atoms:
        new_line = new_line + atom + ' '
    new_line = new_line + '\n'
    final_file.write(new_line)

    # stoichiometry
    new_line = ''
    total_atoms = line_pre.split()
    line_pre = pre_file.readline()
    stoichiometry_atoms = line_pre.split()
    total_num_atoms_stoi = 0
    total_num_atoms = 0
    counter_element = 0
    for atom in diff_atoms:
        total_number_atom = total_atoms.count(atom)
        new_line = new_line + str(int(total_number_atom)*int(stoichiometry_atoms[counter_element])) + ' '
        total_num_atoms_stoi = total_num_atoms_stoi + int(stoichiometry_atoms[counter_element])
        total_num_atoms = total_num_atoms + int(total_number_atom)*int(stoichiometry_atoms[counter_element])
        counter_element = counter_element + 1
    new_line = new_line + '\n'
    final_file.write(new_line)

    atoms_in_poscar = 0
    for element in stoichiometry_atoms:
        atoms_in_poscar = atoms_in_poscar + int(element)

    stoichiometry_unit_cell = []
    for y in range(num_diff_atoms):
        stoichiometry_unit_cell.append(stoichiometry_atoms[y])

    # type of positions
    line_pre = pre_file.readline()
    final_file.write(line_pre)

    pre_file.close()

    # positions
    num_iteration = 1
    for actual_stoi in stoichiometry_unit_cell:
        pre_file = open(f'POSCAR-{x+1:04}_pre', "r")
        for num_lines in range(8):
            pre_file.readline()

        if num_iteration == 1:
            previous_stoi = 0
            next_stoi = total_num_atoms_stoi - int(actual_stoi)
        elif num_iteration == len(stoichiometry_unit_cell):
            previous_stoi = total_num_atoms_stoi - int(actual_stoi)
            next_stoi = 0
        else:
            previous_stoi = 0
            for y in range(num_iteration - 1):
                previous_stoi = previous_stoi + int(stoichiometry_unit_cell[y])
            next_stoi = total_num_atoms_stoi - int(actual_stoi) - previous_stoi
        
        atoms_considered = 0
        while atoms_considered <= atoms_in_poscar:
            if previous_stoi != 0:
                for y in range(previous_stoi):
                    pre_file.readline()
                    atoms_considered = atoms_considered + 1

            for y in range(int(actual_stoi)):
                atom_position_line = pre_file.readline()
                final_file.write(atom_position_line)
                atoms_considered = atoms_considered + 1

            if next_stoi != 0:
                for y in range(next_stoi):
                    pre_file.readline()
                    atoms_considered = atoms_considered + 1

        num_iteration = num_iteration + 1

    # close files
    pre_file.close()
    final_file.close()

# delete the wrong POSCAR files
for x in range(n_structures):
    os.remove(f'POSCAR-{x+1:04}_pre')
"""
