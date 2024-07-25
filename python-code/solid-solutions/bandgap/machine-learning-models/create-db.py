# Pol Benítez Colominas, July 2024
# Universitat Politècnica de Catalunya

# Script to generate a database from DFT calculations

import csv

from pymatgen.io.vasp import Poscar

# Function to read band gaps from DOSCAR file
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

# Save concentrations and bandgaps in arrays
S_conc = []
Br_conc = []
bg = []

path = '../../generate-VCA-grid/VCA_structures/vca-'
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

# Export the data to a csv file
with open('bg-solid-solutions.csv', 'w', newline='') as file:
    writer = csv.writer(file)

    writer.writerow(['S_concentration', 'Br_concentration', 'bandgap'])

    for item1, item2, item3 in zip(S_conc, Br_conc, bg):
        writer.writerow([item1, item2, item3])