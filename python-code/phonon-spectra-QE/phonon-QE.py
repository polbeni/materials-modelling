# Pol Benítez Colominas, October 2024
# University of Cambridge and Universitat Politècnica de Catalunya

# Plot the phonons from the .freq output file provided by QuantumESPRESSO and matdyn.x routine

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

def get_phonon(path):
    """
    Return the phonon dispersion band with structure [[k-point 1], [k-point 2], [k-point 3], ...], where [k-point i] contains all the 
    phonon bands and they are ordered in the decided path in the reciproc space

    Inputs:
        path -> path to de .freq file
    """
    bands = []

    freq_file = open(path, 'r')

    first_line = freq_file.readline()
    num_bands = int(first_line.split()[2].replace(',', ''))
    num_k_points = int(first_line.split()[4])

    for nk in range(num_k_points):
        k_point = []

        freq_file.readline()
        line = freq_file.readline()

        ind_file = 0
        for nb in range(num_bands):
            if (((nb) % 6) == 0) and (nb != 0): 
                line = freq_file.readline()

                ind_file = 0

            k_point.append(float(line.split()[ind_file]))

            ind_file = ind_file + 1

        bands.append(k_point)

    return bands, num_bands, num_k_points

# Get the bands
bands, n_bands, n_kpoints = get_phonon('TiS2.freq')

# Give the proper format to bands, each element is a band, and each band has the different k-points
correct_bands = []
for nb in range(n_bands):
    band = []

    for nk in range(n_kpoints):
        band.append(bands[nk][nb])

    correct_bands.append(band)

# Create an array with the number of k-points
array_k = list(range(0, n_kpoints))

# Interpolate the results
new_array_k = np.linspace(0, n_kpoints, 100)

interp_bands = []
for nb in range(n_bands):
    band = interp1d(array_k, correct_bands[nb], kind='quadratic', fill_value='extrapolate')(new_array_k)

    interp_bands.append(band)

# Plot the results
fig, ax = plt.subplots(figsize=(4, 3))

ax.set_ylim(-150, 450)
ax.axhline(0, linestyle='--', color='black')

ax.set_ylabel('Freq (cm$^{-1}$)')
for nb in range(n_bands):
    ax.plot(new_array_k, interp_bands[nb], color='black')

plt.tight_layout()
plt.savefig('bands.pdf')