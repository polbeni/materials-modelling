# Pol Benítez Colominas, 21/08/2023

import yaml
import numpy as np
import matplotlib.pyplot as plt

def phonons_from_phonopy(file_path):
    """
    This function uses the band.yaml file generated in phonopy calculations to obtain the 
    phonon spectra data. This can be ploted with matplotlib

    Inputs:
        file_path: the path and name of the file

    Outputs:
        num_atoms: the number of atoms in the unit cell
        nqpoints: number of points in the reciprocal space
        npaths: the number of different paths in the reciprocal space
        segments_nqpoints: the number of points in the reciprocal space for each path
        points_labels: the labels of the high symmetry points
        phonons: a matrix array with the information about the phonons, it has the following structure
                    q-point   | branch 1 | branch 2 | branch 3 | ...
                    ------------------------------------------------
                    q-point 1 |   freq   |   freq   |   freq   | ...
                    q-point 2 |   freq   |   freq   |   freq   | ...
                    q-point 3 |   freq   |   freq   |   freq   | ...
                    ...       |   ...    |   ...    |   ...    | ...
    """
    with open(file_path, "r") as file:
        data = yaml.safe_load(file)

    num_atoms = data["natom"]
    nqpoints = data["nqpoint"]
    npaths = data["npath"]
    segments_nqpoints = data["segment_nqpoint"]
    points_labels = data["labels"]

    phonons = np.zeros((nqpoints, num_atoms*3 + 1))
    for x in range(nqpoints):
        phonons[x,0] = data["phonon"][x]["distance"]
        for y in range(num_atoms*3):
            phonons[x,y+1] = data["phonon"][x]["band"][y]["frequency"]

    return num_atoms, nqpoints, npaths, segments_nqpoints, points_labels, phonons

num_atoms, nqpoints, npaths, segments_nqpoints, points_labels, phonons = phonons_from_phonopy('band.yaml')

vertical_lines = ['False']*(npaths-1)
for x in range(npaths-1):
    if points_labels[x][1] != points_labels[x+1][0]:
        vertical_lines[x] = True
    
x_labels = ['point']*(npaths+1)
for x in range(npaths):
    if x == 0:
        x_labels[0] = points_labels[0][0]
    if  x+1 < npaths:
        if vertical_lines[x] == True:
            x_labels[x+1] = points_labels[x][1] + '|' + points_labels[x+1][0]
        else: 
            x_labels[x+1] = points_labels[x][1]
    else:
        x_labels[x+1] = points_labels[x][1]


plt.figure()
plt.xlabel('Reciprocal space')
plt.ylabel('Frequency (THz)')
plt.xlim((phonons[0,0],phonons[-1,0]))
plt.axhline(0, color='black', linestyle='--', linewidth=1)

for x in range(num_atoms*3):
    plt.plot(phonons[:,0], phonons[:,x+1], color='grey', linewidth=.9)

k_point_number = 0
num_segment = 0
distance = 0
for x in segments_nqpoints:
    k_point_number = k_point_number + x
    if vertical_lines[num_segment] == True:
        distance = phonons[k_point_number-1,0]
        plt.axvline(distance, color='black', linewidth=1)
    num_segment = num_segment + 1
    if num_segment >= npaths-1:
        break

x_ticks = [0]*(npaths+1)
k_point_number = 0
element_ticks = 1
for x in segments_nqpoints:
    k_point_number = k_point_number + x
    x_ticks[element_ticks] = phonons[k_point_number-1,0]
    element_ticks = element_ticks + 1

plt.xticks(ticks=x_ticks, labels=x_labels)
plt.savefig('test.pdf')
