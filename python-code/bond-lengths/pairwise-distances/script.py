# Pol Benítez Colominas, October 2025
# The University of Tokyo and Universitat Politècnica de Catalunya

# Compute pairwise distances between all unique atom pairs in a given structure

from pymatgen.core import Structure
from itertools import combinations
from collections import defaultdict

def compute_pairwise_distances(poscar_path, cutoff):
    """
    Compute pairwise distances between all unique atom pairs in a structure

    Inputs:
        poscar_path: path to the POSCAR file
        cutoff: maximum distance to consider for pairwise interactions
    """
    
    structure = Structure.from_file(poscar_path)

    species = [str(sp) for sp in structure.species]

    pair_distances = defaultdict(list)

    # Loop over all unique atom pairs (i < j)
    for i, j in combinations(range(len(structure)), 2):
        sp_i, sp_j = species[i], species[j]

        # Sort species alphabetically for consistent labeling (e.g. 'Ag-S' not 'S-Ag')
        pair_key = "-".join(sorted([sp_i, sp_j]))

        # Compute distance with periodic boundary conditions
        dist = structure.get_distance(i, j)

        if dist <= cutoff:
            pair_distances[pair_key].append(dist)

    # Sort the distance lists
    for key in pair_distances:
        pair_distances[key] = sorted(pair_distances[key])

    return dict(pair_distances)


# Phononic data
path = '/Users/pol/work/upc/ml-gap-prediction/database-CAP/database-uniform/results-PBEsol/pure-compounds/struc-'

Ag_Ag = []
Ag_S = []
Ag_Br = []
S_S = []
Br_S = []
Br_Br = []

for num in range(400):
    struc_path = path + str(num + 1).zfill(4) + '/POSCAR'
    distances = compute_pairwise_distances(struc_path, cutoff=6.0)

    for x in range(len(distances['Ag-Ag'])):
        Ag_Ag.append(distances['Ag-Ag'][x])
    for x in range(len(distances['Ag-S'])):
        Ag_S.append(distances['Ag-S'][x])
    for x in range(len(distances['Ag-Br'])):
        Ag_Br.append(distances['Ag-Br'][x])
    for x in range(len(distances['S-S'])):
        S_S.append(distances['S-S'][x])
    for x in range(len(distances['Br-S'])):
        Br_S.append(distances['Br-S'][x])
    for x in range(len(distances['Br-Br'])):
        Br_Br.append(distances['Br-Br'][x])

key_dataset = 'uniform'

results_uniform = open(f'{key_dataset}_Ag-Ag.txt', 'w')
for x in range(len(Ag_Ag)):
        results_uniform.write(f'{Ag_Ag[x]}\n')
results_uniform.close()

results_uniform = open(f'{key_dataset}_Ag-S.txt', 'w')
for x in range(len(Ag_S)):
        results_uniform.write(f'{Ag_S[x]}\n')
results_uniform.close()

results_uniform = open(f'{key_dataset}_Ag-Br.txt', 'w')
for x in range(len(Ag_Br)):
        results_uniform.write(f'{Ag_Br[x]}\n')
results_uniform.close()

results_uniform = open(f'{key_dataset}_S-S.txt', 'w')
for x in range(len(S_S)):
        results_uniform.write(f'{S_S[x]}\n')
results_uniform.close()

results_uniform = open(f'{key_dataset}_Br-S.txt', 'w')
for x in range(len(Br_S)):
        results_uniform.write(f'{Br_S[x]}\n')
results_uniform.close()

results_uniform = open(f'{key_dataset}_Br-Br.txt', 'w')
for x in range(len(Br_Br)):
        results_uniform.write(f'{Br_Br[x]}\n')
results_uniform.close()


# Uniform data
path = '/Users/pol/work/upc/ml-gap-prediction/database-CAP/database-phononic/results/PBEsol/pure-compounds/struc-'

Ag_Ag = []
Ag_S = []
Ag_Br = []
S_S = []
Br_S = []
Br_Br = []

for num in range(400):
    struc_path = path + str(num + 1).zfill(4) + '/POSCAR'
    distances = compute_pairwise_distances(struc_path, cutoff=6.0)

    for x in range(len(distances['Ag-Ag'])):
        Ag_Ag.append(distances['Ag-Ag'][x])
    for x in range(len(distances['Ag-S'])):
        Ag_S.append(distances['Ag-S'][x])
    for x in range(len(distances['Ag-Br'])):
        Ag_Br.append(distances['Ag-Br'][x])
    for x in range(len(distances['S-S'])):
        S_S.append(distances['S-S'][x])
    for x in range(len(distances['Br-S'])):
        Br_S.append(distances['Br-S'][x])
    for x in range(len(distances['Br-Br'])):
        Br_Br.append(distances['Br-Br'][x])

key_dataset = 'phononic'

results_uniform = open(f'{key_dataset}_Ag-Ag.txt', 'w')
for x in range(len(Ag_Ag)):
        results_uniform.write(f'{Ag_Ag[x]}\n')
results_uniform.close()

results_uniform = open(f'{key_dataset}_Ag-S.txt', 'w')
for x in range(len(Ag_S)):
        results_uniform.write(f'{Ag_S[x]}\n')
results_uniform.close()

results_uniform = open(f'{key_dataset}_Ag-Br.txt', 'w')
for x in range(len(Ag_Br)):
        results_uniform.write(f'{Ag_Br[x]}\n')
results_uniform.close()

results_uniform = open(f'{key_dataset}_S-S.txt', 'w')
for x in range(len(S_S)):
        results_uniform.write(f'{S_S[x]}\n')
results_uniform.close()

results_uniform = open(f'{key_dataset}_Br-S.txt', 'w')
for x in range(len(Br_S)):
        results_uniform.write(f'{Br_S[x]}\n')
results_uniform.close()

results_uniform = open(f'{key_dataset}_Br-Br.txt', 'w')
for x in range(len(Br_Br)):
        results_uniform.write(f'{Br_Br[x]}\n')
results_uniform.close()