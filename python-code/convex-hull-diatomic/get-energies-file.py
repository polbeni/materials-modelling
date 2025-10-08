# Pol Benítez Colominas, October 2025
# The University of Tokyo and Universitat Politècnica de Catalunya

# Generates the files with energies from CrySPY calculations


import csv


def get_formula(pair):
    """
    Returns the formula for a given pair ready for pymatgen phase diagram function
    (IMPORTANT: limited to a pairs below 10, with 10 not included)

    Inputs:
        pair: stoichiometry of AB compund
    """

    formula = ''

    stoi_A = pair[0]
    stoi_B = pair[2]

    formula = formula + 'A'
    if stoi_A != '1':
        formula = formula + stoi_A
    
    formula = formula + 'B'
    if stoi_B != '1':
        formula = formula + stoi_B

    return formula


# Get the pairs information
pairs = []

with open('pairs.txt', 'r') as file:
    for line in file:
        pairs.append(line.split()[0])

# Get the formulas in format for pymatgen
formulas = []

for pair in pairs:
    formulas.append(get_formula(pair))

# Get the energies (above a threshold)
energy_threshold = -4.0 # eV/atom, below this it's likely to be an error by the MLIPs

calculations_results = 'calculations/'

data_structure = []

for pair in range(len(pairs)):
    pair_path = calculations_results + pairs[pair] + '/data/cryspy_rslt_energy_asc'

    with open(pair_path, 'r') as file:
        file.readline()

        for line in file:
            if float(line.split()[6]) > energy_threshold:
                data_structure.append([formulas[pair], float(line.split()[6])])

# Manually add results for the primary compounds
data_structure.append(['A', -1])
data_structure.append(['B', -1])

# Save it in a csv file
with open('energies.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    
    writer.writerow(['formula', 'energy_per_atom'])
    
    writer.writerows(data_structure)
