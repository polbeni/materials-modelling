# Pol Benítez Colominas, October 2024
# University of Cambridge and Universitat Politècnica de Catalunya

# Applies gaussian noise to distort a given structure

import numpy as np

from pymatgen.io.vasp import Poscar


#### Basic variables definition ####
noise = 'gaussian'          # types of noise: gaussian or uniform

width_disp = 0.5            # standard deviation (width) of the gaussian noise
max_disp = 0.1              # maximum displacement for uniform noise

num_gen_struc = 10          # number of generated structures

path_to_save = 'results'    # path to save the resulting structures


#### Distorted structures generation ####
poscar = Poscar.from_file("POSCAR-2x2x2")
original_poscar = poscar.structure

for num_struc in range(num_gen_struc):
    cp_poscar = original_poscar.copy()

    for site in cp_poscar:
        if noise == 'gaussian':
            displacement_vector = np.random.uniform(0, width_disp, 3)
        elif noise == 'uniform':
            displacement_vector = np.random.normal(-max_disp, max_disp, 3)

        site.coords = site.coords + displacement_vector

    path_struc = path_to_save + '/POSCAR-' + str(num_struc + 1).zfill(5)
    cp_poscar.to(filename=path_struc, fmt='Poscar')
