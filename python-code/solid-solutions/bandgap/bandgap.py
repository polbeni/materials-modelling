import numpy as np
import matplotlib.pyplot as plt

from pymatgen.io.vasp import Poscar

def electronic_bandGap(file_name):
    """
    This functions uses DOSCAR file generated in VASP simulations and returns the Fermi energy
    the band gap, and the energies of the band gap (respect the exchange-correlation functional
    used).

    file_name: path of the DOSCAR file
    """
    
    file = open(file_name, "r")

    for x in range(6):
        actual_string = file.readline()
        if x == 5:
            fermiEnergy = float(actual_string.split()[3])

    file.close()

    file = open(file_name, "r")

    for x in range(6):
        file.readline()

    for x in file:
        actual_string = x

        if (float(actual_string.split()[0]) <= fermiEnergy+0.1) and (float(actual_string.split()[0]) >= fermiEnergy-0.1):
            density_bandGap = float(actual_string.split()[2])

            break

    file.close()

    file = open(file_name, "r")

    for x in range(6):
        file.readline()

    for x in file:
        actual_string = x

        if float(actual_string.split()[2]) == density_bandGap:
            minEnergy = float(actual_string.split()[0])

            break   

    for x in file:
        actual_string = x

        if float(actual_string.split()[2]) != density_bandGap:
            maxEnergy = float(actual_string.split()[0])

            break 
    bandGap = maxEnergy - minEnergy

    file.close()
    
    return fermiEnergy, minEnergy, maxEnergy, bandGap

S_conc = []
Br_conc = []
bg = []

path = '../generate-VCA-grid/VCA_structures/vca-'
for x in range(121):
    path_poscar = path + str(x + 1).zfill(3) + '/bandgap/POSCAR'
    poscar = Poscar.from_file(path_poscar)
    structure = poscar.structure
    list_atoms = [site.species_string for site in structure.sites]
    unique_atom_list = list(dict.fromkeys(list_atoms))

    path_incar = path + str(x + 1).zfill(3) + '/bandgap/INCAR'
    coefs = []
    with open(path_incar, 'r') as file:
        for line in file:
            values = line.split()
    for y in range(len(values)):
        if (y != 0) and (y != 1):
            coefs.append(float(values[y]))
    if 'S' in unique_atom_list:
        S_conc.append(coefs[1])
    else:
        S_conc.append(0)
    if 'Br' in unique_atom_list:
        if ('S' in unique_atom_list) and ('Se' in unique_atom_list):
            Br_conc.append(coefs[3])
        if ('S' not in unique_atom_list) or ('Se' not in unique_atom_list):
            Br_conc.append(coefs[2])
    else:
        Br_conc.append(0)

    path_bg = path + str(x + 1).zfill(3) + '/bandgap/DOSCAR'
    _, _, _, bandgap = electronic_bandGap(path_bg)
    bg.append(bandgap)
    

matrix_bg = []
iteration = 0
for x in range(11):
    row = []
    for y in range(11):
        row.append(bg[120 - iteration])

        iteration = iteration + 1
    
    n_row = []
    for y in range(11):
        n_row.append(row[10 - y])
    
    matrix_bg.append(n_row)

fig, ax = plt.subplots(figsize=(6,4))

cax = ax.imshow(matrix_bg, cmap='magma', aspect='auto')
cbar = fig.colorbar(cax, ax=ax)
cbar.set_label('$E_{g}$ (eV)', fontsize=12) 

ax.set_title('Ag$_3$S$_{x}$Se$_{1-x}$Br$_{y}$I$_{1-y}$')
ax.set_xlabel('x', fontsize=12)
ax.set_ylabel('y', fontsize=12)

x_ticks = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
x_labels = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 10]
ax.set_xticks(ticks=x_ticks, labels=x_labels)

y_ticks = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
y_labels = [1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0]
ax.set_yticks(ticks=y_ticks, labels=y_labels)

plt.tight_layout()
plt.savefig('solid-solutions-bg.pdf')