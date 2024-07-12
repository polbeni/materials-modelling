# Pol Benítez Colominas, July 2024
# Universitat Politècnica de Catalunya

# Generate a distance heat map of atoms distance after a distortion

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

from pymatgen.core.structure import Structure

def length_bond(struc, atom1, atom2):
    """
    Computes the lenght of the bond between two atoms

    Inputs:
        struc -> unit cell structure in pymatgen format
        atom1 -> first atom index
        atom2 -> second atom index
    """

    length = struc.get_distance(atom1, atom2)

    return length

def change_distance(d_orginal, d_distorted):
    """
    Computes the length variation after a distortion as (d-d')/d (where d is orginal distance and d' distorted)
    If the variation >0 the two atoms are closer, while <0 if atoms are further

    Inputs:
        d_original: original distance
        d_distorted: distorted distance
    """

    variation = ((d_orginal - d_distorted)/d_orginal)*100

    return variation

#################### Ag3SBr ####################

names_bond = ['S-Ag$_1$', 'S-Ag$_2$', 'S-Ag$_3$', 'S-Br', 'Br-Ag$_1$', 'Br-Ag$_2$', 'Br-Ag$_3$', 
              'Ag$_1$-Ag$_2$', 'Ag$_1$-Ag$_3$', 'Ag$_2$-Ag$_3$']

original_distance = [] # 0:S-Ag1, 1:S-Ag2, 2:S-Ag3, 3:S-Br,  4:Br-Ag1, 5:Br-Ag2, 6:Br-Ag3, 7:Ag1-Ag2, 8:Ag1-Ag3, 9:Ag2-Ag3

poscar_path = '/home/pol/work/materials-modelling/materials-simulations/chalcohalide-antiperovskites/silver/Ag-Pm3m/bg-phonon-effect/electronic-structure-Pm3m/POSCAR'
structure_original = Structure.from_file(poscar_path)

original_distance.append(length_bond(structure_original, 3, 0))
original_distance.append(length_bond(structure_original, 3, 1))
original_distance.append(length_bond(structure_original, 3, 2))
original_distance.append(length_bond(structure_original, 3, 4))
original_distance.append(length_bond(structure_original, 4, 0))
original_distance.append(length_bond(structure_original, 4, 1))
original_distance.append(length_bond(structure_original, 4, 2))
original_distance.append(length_bond(structure_original, 0, 1))
original_distance.append(length_bond(structure_original, 0, 2))
original_distance.append(length_bond(structure_original, 1, 2))

var_phonon = []
path_struct = '/home/pol/work/materials-modelling/materials-simulations/chalcohalide-antiperovskites/silver/Ag-Pm3m/bg-phonon-effect/distorted-structures/'
for x in range(15):
    phonon_dist = []

    phonon_poscar_path = path_struct + 'phonon' + str(x + 1).zfill(2) + '-4/POSCAR'
    structure_phonon = Structure.from_file(phonon_poscar_path)

    phonon_dist.append(length_bond(structure_phonon, 3, 0))
    phonon_dist.append(length_bond(structure_phonon, 3, 1))
    phonon_dist.append(length_bond(structure_phonon, 3, 2))
    phonon_dist.append(length_bond(structure_phonon, 3, 4))
    phonon_dist.append(length_bond(structure_phonon, 4, 0))
    phonon_dist.append(length_bond(structure_phonon, 4, 1))
    phonon_dist.append(length_bond(structure_phonon, 4, 2))
    phonon_dist.append(length_bond(structure_phonon, 0, 1))
    phonon_dist.append(length_bond(structure_phonon, 0, 2))
    phonon_dist.append(length_bond(structure_phonon, 1, 2))

    phonon_var = []
    for y in range(10):
        phonon_var.append(change_distance(original_distance[y], phonon_dist[y]))
    
    var_phonon.append(phonon_var)

var_phonon_np = np.array(var_phonon)
var_phonon_np = var_phonon_np.T


fig, ax = plt.subplots(figsize=(6,4))

cax = ax.imshow(var_phonon_np, cmap='magma', aspect='auto')
cbar = fig.colorbar(cax, ax=ax)
cbar.set_label('Bond variation (%)')  # Set your label here

ax.invert_yaxis()

ax.set_title('Ag$_3$SBr')
ax.set_xlabel('Phonon mode #')
ax.set_ylabel('Bond')

x_ticks = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]
x_labels = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
ax.set_xticks(ticks=x_ticks, labels=x_labels)

y_ticks = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
y_labels = names_bond
ax.set_yticks(ticks=y_ticks, labels=y_labels)

major_locator = MultipleLocator(1)  
ax.xaxis.set_major_locator(major_locator)

plt.tight_layout()
plt.savefig('bond_variation-Ag3SBr.pdf')
################################################


#################### BaTiO3 ####################

names_bond = ['Ti-O$_1$', 'Ti-O$_2$', 'Ti-O$_3$', 'Ti-Ba', 'Ba-O$_1$', 'Ba-O$_2$', 'Bo-O$_3$', 
              'O$_1$-O$_2$', 'O$_1$-O$_3$', 'O$_2$-O$_3$']

original_distance = [] # 0:Ti-O1, 1:Ti-O2, 2:Ti-O3, 3:Ti-Ba,  4:Ba-O1, 5:Ba-O2, 6:Ba-O3, 7:O1-O2, 8:O1-O3, 9:O2-O3

poscar_path = '/home/pol/work/materials-modelling/materials-simulations/perovskites-ep-coupling/simulations-converged/mp-2998-BaTiO3/electronic-structure/dist-00/POSCAR'
structure_original = Structure.from_file(poscar_path)

original_distance.append(length_bond(structure_original, 3, 0))
original_distance.append(length_bond(structure_original, 3, 1))
original_distance.append(length_bond(structure_original, 3, 2))
original_distance.append(length_bond(structure_original, 3, 4))
original_distance.append(length_bond(structure_original, 4, 0))
original_distance.append(length_bond(structure_original, 4, 1))
original_distance.append(length_bond(structure_original, 4, 2))
original_distance.append(length_bond(structure_original, 0, 1))
original_distance.append(length_bond(structure_original, 0, 2))
original_distance.append(length_bond(structure_original, 1, 2))

var_phonon = []
path_struct = '/home/pol/work/materials-modelling/materials-simulations/perovskites-ep-coupling/BaTiO3-phonon-effect/bg-phonon-effect/'
for x in range(15):
    phonon_dist = []

    phonon_poscar_path = path_struct + 'phonon' + str(x + 1).zfill(2) + '-4/POSCAR'
    structure_phonon = Structure.from_file(phonon_poscar_path)

    phonon_dist.append(length_bond(structure_phonon, 3, 0))
    phonon_dist.append(length_bond(structure_phonon, 3, 1))
    phonon_dist.append(length_bond(structure_phonon, 3, 2))
    phonon_dist.append(length_bond(structure_phonon, 3, 4))
    phonon_dist.append(length_bond(structure_phonon, 4, 0))
    phonon_dist.append(length_bond(structure_phonon, 4, 1))
    phonon_dist.append(length_bond(structure_phonon, 4, 2))
    phonon_dist.append(length_bond(structure_phonon, 0, 1))
    phonon_dist.append(length_bond(structure_phonon, 0, 2))
    phonon_dist.append(length_bond(structure_phonon, 1, 2))

    phonon_var = []
    for y in range(10):
        phonon_var.append(change_distance(original_distance[y], phonon_dist[y]))
    
    var_phonon.append(phonon_var)

var_phonon_np = np.array(var_phonon)
var_phonon_np = var_phonon_np.T


fig, ax = plt.subplots(figsize=(6,4))

cax = ax.imshow(var_phonon_np, cmap='magma', aspect='auto')
cbar = fig.colorbar(cax, ax=ax)
cbar.set_label('Bond variation (%)')  # Set your label here

ax.invert_yaxis()

ax.set_title('BaTiO$_3$')
ax.set_xlabel('Phonon mode #')
ax.set_ylabel('Bond')

x_ticks = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]
x_labels = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
ax.set_xticks(ticks=x_ticks, labels=x_labels)

y_ticks = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
y_labels = names_bond
ax.set_yticks(ticks=y_ticks, labels=y_labels)

major_locator = MultipleLocator(1)  
ax.xaxis.set_major_locator(major_locator)

plt.tight_layout()
plt.savefig('bond_variation-BaTiO3.pdf')
################################################