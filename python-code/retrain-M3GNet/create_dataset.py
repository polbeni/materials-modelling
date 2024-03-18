# Pol Benítez Colominas, March 2024
# Universitat Politècnica de Catalunya

# Script to create a dataset with DFT data to re-train M3GNet
# IMPORTANT: this code is based in the code developed by Cibrán: https://github.com/CibranLopez/m3gnet

import os
import warnings

import pandas as pd

from pymatgen.io.vasp.outputs import Vasprun

warnings.simplefilter('ignore')

# define the materials with DFT data, in this I have a directory of data, that conatains folders for each
# of the materials, and each materials contains a different number of folders with vasprun.xml files,
# for example: data_training/Ag3SCl/01/vasprun.xml
materials = ['Ag3SCl', 'Ag3SBr', 'Ag3SI', 'Ag3SeCl', 'Ag3SeBr', 'Ag3SeI', 
             'Cu3SCl', 'Cu3SBr', 'Cu3SI', 'Cu3SeCl', 'Cu3SeI']

data_path = 'data_training/'

# we want to save four features for each ionic step, the structure, the energy, the forces and the stresses
structures = []
energies = []
forces = []
stresses = []

# generate a file where we store the number of data rows (accumulated) for each material
info_data = open('info_data.txt', 'w')
info_data.write('Material   #rows\n')
num_rows = 0

# run for each material
for mat in materials:
    phases = [d for d in os.listdir(data_path + mat) if os.path.isdir(os.path.join(data_path + mat, d))]

    # run for each phase of the given material
    for phase in phases:
        try:
            vasprun = Vasprun(data_path + mat + '/' + phase + '/vasprun.xml', exception_on_bad_xml=False)
        except:
            print('Error: vasprun not correctly loaded.')
            continue

        # save the desired features for each ionic step of each phase of the given material
        for step in vasprun.ionic_steps:
            structures.append(step['structure'])
            energies.append(step['electronic_steps'][-1]['e_fr_energy'])
            forces.append(step['forces'])
            stresses.append(step['stress'])

            num_rows = num_rows + 1
    
    info_data.write(f'{mat}   {num_rows:06d}\n')

info_data.close()

# save the date in a pandas dataframe object and save it in a csv file
data = {
    'structure': structures,
    'energy': energies,
    'force': forces,
    'stress': stresses
}

df_data = pd.DataFrame(data)
df_data.to_csv('data_relaxation.csv', index=False)