# Pol Benítez Colominas, February 2025
# Universitat Politècnica de Catalunya

# Generates the pT diagram from the phonopy-qha results

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

def get_gibbs_energy(file_path, num_atoms):
    """
    This function extracts the temperature and gibbs energy (per atom) from the gibbs energy of phonopy-qha

    Inputs:
        file_path: the path and name of the file
        num_atoms: num atoms in the unit cell
    """

    file = open(file_path, 'r')

    T = []
    G = []

    for x in range(181):
        line = file.readline()
        T.append(float(line.split()[0]))
        G.append(((float(line.split()[1]) / (num_atoms)))*3)

    return T, G


############ ANHARMONIC ############

phases = ['monoGS', 'cubic', 'orthoAPref', 'orthoI', 'orthoIIIFE', 'orthoIstar', 'tetra']
phases_dict = {
    'monoGS': 1,
    'orthoIstar': 2,
    'orthoI': 3,
    'orthoIIIFE': 4,
    'orthoAPref': 5,
    'tetra': 6,
    'cubic': 7
}
num_atoms_phases = [12, 12, 12, 24, 12, 24, 12]

pressure_points = 200
temperature_points = 181

pressure_array = np.linspace(0, 10, pressure_points) # in GPa
temperature_array = np.linspace(0, 1800, temperature_points) # in K

less_energy = np.zeros([pressure_points, temperature_points])

for pressure in range(pressure_points):
    for temperature in range(temperature_points):
        _, gibbs = get_gibbs_energy(f'data-pressure/monoGS/gibbs-{pressure + 1}.dat', 12)
        energy_monoGS = gibbs[temperature]

        _, gibbs = get_gibbs_energy(f'data-pressure/orthoIstar/gibbs-{pressure + 1}.dat', 24)
        energy_orthoIstar = gibbs[temperature]

        _, gibbs = get_gibbs_energy(f'data-pressure/orthoI/gibbs-{pressure + 1}.dat', 24)
        energy_orthoI = gibbs[temperature]

        _, gibbs = get_gibbs_energy(f'data-pressure/orthoIIIFE/gibbs-{pressure + 1}.dat', 12)
        energy_orthoIIIFE = gibbs[temperature]

        _, gibbs = get_gibbs_energy(f'data-pressure/orthoAPref/gibbs-{pressure + 1}.dat', 12)
        energy_orthoAPref = gibbs[temperature]

        _, gibbs = get_gibbs_energy(f'data-pressure/tetra/gibbs-{pressure + 1}.dat', 12)
        energy_tetra = gibbs[temperature]

        _, gibbs = get_gibbs_energy(f'data-pressure/cubic/gibbs-{pressure + 1}.dat', 12)
        energy_cubic = gibbs[temperature]

        smallest_value = min(energy_monoGS, energy_orthoIstar, energy_orthoI, energy_orthoIIIFE, energy_orthoAPref, energy_tetra, energy_cubic)

        if energy_monoGS == smallest_value:
            less_energy[pressure, temperature] = phases_dict['monoGS']
        elif energy_orthoIstar == smallest_value:
            less_energy[pressure, temperature] = phases_dict['orthoIstar']
        elif energy_orthoI == smallest_value:
            less_energy[pressure, temperature] = phases_dict['orthoI']
        elif energy_orthoIIIFE == smallest_value:
            less_energy[pressure, temperature] = phases_dict['orthoIIIFE']
        elif energy_orthoAPref == smallest_value:
            less_energy[pressure, temperature] = phases_dict['orthoAPref']
        elif energy_tetra == smallest_value:
            less_energy[pressure, temperature] = phases_dict['tetra']
        elif energy_cubic == smallest_value:
            less_energy[pressure, temperature] = phases_dict['cubic']

less_energy = less_energy[::-1, :]

fig, ax = plt.subplots(figsize=(4, 3))

ax.set_title('HfO$_2$ phase diagram') 

ax.set_xlabel('Temperature (K)')
ax.set_ylabel('Pressure (GPa)')

im = ax.imshow(less_energy, cmap='plasma', aspect='auto')

x_ticks = [0, 90, 179]
x_labels = [0, 900, 1800]

y_ticks = [0, 99, 199] 
y_labels = [10, 5, 0]  

ax.set_xticks(x_ticks, x_labels)
ax.set_yticks(y_ticks, y_labels)

ax.text(70, 150, 'monoGS', color='white')
ax.text(70, 50, 'orthoI', color='black')
ax.text(150, 50, 'tetra', color='black')

plt.tight_layout()
plt.savefig('pt-diagram.pdf')
