# Pol Benítez Colominas, June 2025
# Universitat Politècnica de Catalunya

# Group all the phases of the same material (chemical composition) and group by number of chemical species

import os
import shutil

# Read and save all the data
materials = []
chemical_structures = []

with open('all_structures/materials_info.txt', 'r') as file:
    file.readline()

    lines = file.readlines()
    for line in lines:
        materials.append([line.split()[0], line.split()[1].replace('(', '_').replace(')', '_'), line.split()[2]])
        chemical_structures.append(line.split()[1].replace('(', '_').replace(')', '_'))


# Create a set with non repeated chemical structures
unique_structures = set(chemical_structures)
print('The total number of structures is:            ', len(chemical_structures))
print('The total number of unique materials is:      ', len(unique_structures))

# Group all the materials with more than one phase
grouped_compounds = []
num_unique = 1
for material in unique_structures:
    grouped = []
    for structure in range(len(materials)):
        if material == materials[structure][1]:
            num_chemical = materials[structure][2]

            grouped.append(materials[structure][0])
    
    if len(grouped) > 1:
        grouped_compounds.append([material, num_chemical, grouped])

    print(f'{num_unique} structures checked of a total of {len(unique_structures)}')
    print(f'{len(grouped_compounds)} materials with more than one phase')

    num_unique = num_unique + 1

# Copy the structures in different folders
path_dir = 'structures-grouped'
if os.path.exists(path_dir):
    shutil.rmtree(path_dir)
os.mkdir(path_dir)

for num_comp in range(9):
    path_dir = f'structures-grouped/{num_comp + 1}-materials'
    if os.path.exists(path_dir):
        shutil.rmtree(path_dir)
    os.mkdir(path_dir)

    for group in grouped_compounds:
        if int(group[1]) == num_comp + 1:
            path_dir = f'structures-grouped/{num_comp + 1}-materials/{group[0]}'
            if os.path.exists(path_dir):
                shutil.rmtree(path_dir)
            os.mkdir(path_dir)

            for mat in group[2]:
                shutil.copy('all_structures/' + mat + '.cif', path_dir + '/' + mat + '.cif')