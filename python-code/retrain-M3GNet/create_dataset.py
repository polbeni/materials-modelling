# Pol Benítez Colominas, March 2024
# Universitat Politècnica de Catalunya

# Script to create a dataset with DFT data to re-train M3GNet
# IMPORTANT: this code is based in the code developed by Cibrán: https://github.com/CibranLopez/m3gnet

import os
import warnings

import pandas as pd

from pymatgen.io.vasp.outputs import Vasprun

warnings.simplefilter('ignore')


materials = ['Ag3SCl', 'Ag3SBr', 'Ag3SI', 'Ag3SeCl', 'Ag3SeBr', 'Ag3SeI', 
             'Cu3SCl', 'Cu3SBr', 'Cu3SI', 'Cu3SeCl', 'Cu3SeI']

data_path = 'data_training/'

columns = ['structure', 'energy', 'force', 'stress']
structures = []
energies = []
forces = []
stresses = []

for mat in materials:
    phases = [d for d in os.listdir(data_path + mat) if os.path.isdir(os.path.join(data_path + mat, d))]

    for phase in phases:
        try:
            vasprun = Vasprun(data_path + mat + '/' + phase + '/vasprun.xml', exception_on_bad_xml=False)
        except:
            print('Error: vasprun not correctly loaded.')
            continue

        for step in vasprun.ionic_steps:
            structures.append(step['structure'])
            energies.append(step['electronic_steps'][-1]['e_fr_energy'])
            forces.append(step['forces'])
            stresses.append(step['stress'])

data = {
    'structure': structures,
    'energy': energies,
    'force': forces,
    'stress': stresses
}

df_data = pd.DataFrame(data)
df_data.to_csv('data_relaxation.csv', index=False)