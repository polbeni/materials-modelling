# Pol Benítez Colominas, October 2025
# The University of Tokyo and Universitat Politècnica de Catalunya

# Generate calculations for all the possible solid solutions

import os
import shutil
from math import gcd


### FUNCTIONS ###
def unique_pairs(n):
    """
    Generates all the possible distinct pairs (i,j) with no multiplicity with 1<=i,j<=n

    Inputs:
        n: maximum number for the pair generation
    """

    pairs = []
    for i in range(1, n + 1):
        for j in range(1, n + 1):
            if gcd(i, j) == 1:
                pairs.append((i, j))

    return pairs

def file_pairs(pairs, file_name):
    """
    Creates a files with all the possible pairs

    Inputs:
        pairs: list with all the generated pairs
        file_name: name of the generated file
    """

    file = open(file_name, 'w')
    
    for pair in pairs:
        file.write(f'{pair[0]}_{pair[1]}\n')

    file.close()

def generate_cryspy_in_file(path_original, pair, path_to_save):
    """
    Generated the cryspy.in file with the pair number of atoms

    Inputs:
        path_original: path of the original cryspy.in file
        pair: pair of the stoichimetry
        path_to_save: path to save the file
    """

    original_file = open(path_original, 'r')
    new_file = open(path_to_save, 'w')

    for _ in range(10):
        line = original_file.readline()
        new_file.write(line)

    original_file.readline()
    new_file.write(f'nat = {pair[0]} {pair[1]}\n')

    for _ in range(17):
        line = original_file.readline()
        new_file.write(line)

    original_file.close()
    new_file.close()


### MAIN ###
# Generate pairs and save as file
pairs = unique_pairs(8)
file_pairs(pairs, 'pairs.txt')

# Create a dir to do calculations
path_name = 'calculations'

if os.path.exists(path_name):
    shutil.rmtree(path_name)
os.makedirs(path_name)

# Generate paths ready to be submited
bash_file = open(f'{path_name}/run-cryspy.sh', 'w')
bash_file.write('#!/bin/bash\n')
bash_file.write('\n')

times_bash = 200
delay_bash = 5 # in seconds

bash_file.write(f'TIMES={times_bash}\n')
bash_file.write(f'DELAY={delay_bash}\n')
bash_file.write('\n')

for pair in pairs:
    path_pair = path_name + '/' + f'{pair[0]}_{pair[1]}' + '/'
    os.makedirs(path_pair)

    shutil.copytree('cryspy-files/calc_in', path_pair + 'calc_in')

    generate_cryspy_in_file('cryspy-files/cryspy.in', pair, path_pair + 'cryspy.in')

    bash_file.write(f'cd {pair[0]}_{pair[1]}\n')
    bash_file.write('for i in $(seq 1 $TIMES); do\n')
    bash_file.write('    cryspy\n')
    bash_file.write('    sleep $DELAY\n')
    bash_file.write('done\n')
    bash_file.write('cd ..\n')
    bash_file.write('\n')

bash_file.close()

os.system(f'chmod +x {path_name}/run-cryspy.sh')