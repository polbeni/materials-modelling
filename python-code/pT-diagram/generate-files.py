# Pol Benítez Colominas, February 2025
# Universitat Politècnica de Catalunya

# Generates an sh file that will execute phonopy-qha for a range of T and p, and the desired phases

import os
import numpy as np

pressure_array = np.linspace(0, 10, 200) # in GPa

label = open('pressure_lable.txt', 'w')
for pressure in range(len(pressure_array)):
    label.write(f'{pressure + 1}     {pressure_array[pressure]}\n')

sh_script = open('generate-g-pressure.sh', 'w')
sh_script.write('#!/bin/bash \n')
sh_script.write(' \n')

phases_list = ['monoGS', 'orthoIstar', 'orthoI', 'orthoIIIFE', 'orthoAPref', 'tetra', 'cubic']

for phase in phases_list:
    sh_script.write(f'cd data-pressure/{phase} \n')

    sh_script.write(f'cp ../../data/{phase}/e-v.dat . \n')
    sh_script.write(f'cp ../../data/{phase}/thermal_properties-0.yaml . \n')
    sh_script.write(f'cp ../../data/{phase}/thermal_properties-1.yaml . \n')
    sh_script.write(f'cp ../../data/{phase}/thermal_properties-2.yaml . \n')
    sh_script.write(f'cp ../../data/{phase}/thermal_properties-3.yaml . \n')
    sh_script.write(f'cp ../../data/{phase}/thermal_properties-4.yaml . \n')

    for pressure in range(len(pressure_array)):
        sh_script.write(f'phonopy-qha --tmax 1800 --pressure {pressure_array[pressure]} e-v.dat thermal_properties-*.yaml \n')
        sh_script.write(f'cp gibbs-temperature.dat gibbs-{pressure + 1}.dat \n')

    sh_script.write(f'rm Cp-temperature.dat Cp-temperature_polyfit.dat Cv-volume.dat bulk_modulus-temperature.dat dsdv-temperature.dat entropy-volume.dat gibbs-temperature.dat gruneisen-temperature.dat helmholtz-volume.dat helmholtz-volume_fitted.dat thermal_expansion.dat volume-temperature.dat \n')

    sh_script.write('cd ../.. \n')

os.system("chmod +x generate-g-pressure.sh")
