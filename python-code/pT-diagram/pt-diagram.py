# Pol Benítez Colominas, February 2025
# Universitat Politècnica de Catalunya

# Generates the pT diagram from the phonopy-qha results

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib.colors as mcolors

def get_gibbs_energy(file_path, num_atoms):
    """
    This function extracts the temperature and gibbs energy (per atom) from the gibbs energy of phonopy-qha

    Inputs:
        file_path -> the path and name of the file
        num_atoms -> num atoms in the unit cell
    """

    file = open(file_path, 'r')

    T = []
    G = []

    for x in range(181):
        line = file.readline()
        T.append(float(line.split()[0]))
        G.append(((float(line.split()[1]) / (num_atoms)))*3)

    return T, G


def read_volume(path, temperature_in):
    """
    It reads the value volume for a given pressure and temperature from the volume-temperature.dat file

    Inputs:
        path -> path to the volume-temperature.dat file
        temperature_in -> desired temperature (its index)
    """

    vol_temp = open(path, 'r')
    for _ in range(temperature_in + 1):
        line = vol_temp.readline()
    volume = float(line.split()[1])
    vol_temp.close()

    return volume


def check_valid_T(pressure_in, temperature_in):
    """
    It verifies the minimum T at what the given p computations are valid (no extrapolated)
    It looks for all the phases at takes the maximum minimum T of the phases

    Inputs: 
        pressure_in -> desired pressure (its index)
        temperature_in -> desired temperature (its index)
    """
    
    # define condition, True if there is a volume smaller, False otherwise
    condition = False 

    val_monoGS = read_volume(f'data-pressure/monoGS/volume-temperature-{pressure_in + 1}.dat', temperature_in)
    if val_monoGS <= 132.80:
        condition = True

    val_cubic = read_volume(f'data-pressure/cubic/volume-temperature-{pressure_in + 1}.dat', temperature_in)
    if val_cubic <= 123.90:
        condition = True

    val_orthoAPref = read_volume(f'data-pressure/orthoAPref/volume-temperature-{pressure_in + 1}.dat', temperature_in)
    if val_orthoAPref <= 141.68:
        condition = True

    val_orthoI = read_volume(f'data-pressure/orthoI/volume-temperature-{pressure_in + 1}.dat', temperature_in)
    if val_orthoI <= 254.50:
        condition = True

    val_orthoIIIFE = read_volume(f'data-pressure/orthoIIIFE/volume-temperature-{pressure_in + 1}.dat', temperature_in)
    if val_orthoIIIFE <= 128.00:
        condition = True

    val_orthoIstar = read_volume(f'data-pressure/orthoIstar/volume-temperature-{pressure_in + 1}.dat', temperature_in)
    if val_orthoIstar <= 266.20:
        condition = True

    val_tetra = read_volume(f'data-pressure/tetra/volume-temperature-{pressure_in + 1}.dat', temperature_in)
    if val_tetra <= 126.70:
        condition = True

    return condition


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

        extrapolated = check_valid_T(pressure, temperature)
        if extrapolated == True:
            less_energy[pressure, temperature] = 0

less_energy = less_energy[::-1, :]

fig, ax = plt.subplots(figsize=(4, 3))

ax.set_title('HfO$_2$ phase diagram (anharmonic)') 

ax.set_xlabel('Temperature (K)')
ax.set_ylabel('Pressure (GPa)')


colors = ['white', 'bisque', 'green', 'turquoise', 'yellow', 'purple', 'mediumpurple', 'purple']
values = [0, 1, 2, 3, 4, 5, 6, 7] 

cmap = mcolors.ListedColormap(colors)
norm = mcolors.BoundaryNorm(boundaries=[-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5], ncolors=len(values))

im = ax.imshow(less_energy, cmap=cmap, norm=norm, aspect='auto')

x_ticks = [0, 90, 179]
x_labels = [0, 900, 1800]

y_ticks = [0, 99, 199] 
y_labels = [10, 5, 0]  

ax.set_xticks(x_ticks, x_labels)
ax.set_yticks(y_ticks, y_labels)

ax.text(70, 150, 'monoGS', color='black')
ax.text(90, 75, 'orthoI', color='black')
ax.text(150, 60, 'tetra', color='black')
ax.text(30, 40, '?????', color='black')

plt.tight_layout()
plt.savefig('pt-diagram.pdf')
