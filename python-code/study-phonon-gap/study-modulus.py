import numpy as np
import matplotlib.pyplot as plt

def compute_modulus(x, y, z):
    """
    Computes the modulus of a vector

    Inputs:
        x, y, z -> cartesian components of the vector
    """

    return (x**2 + y**2 + z**2)**(1 / 2)

OUTCAR = open('gamma-phonons/OUTCAR', 'r')
for x in range(134096):
    OUTCAR.readline()

tot_mod = []
for _ in range(60): # phonon loop
    OUTCAR.readline()
    OUTCAR.readline()

    phonon = []
    for _ in range(20): # atom loop
        line = OUTCAR.readline()
        phonon.append(compute_modulus(float(line.split()[3]), float(line.split()[4]), float(line.split()[5])))
    
    line = OUTCAR.readline()

    tot_mod.append(phonon)

OUTCAR.close()

matrix_modulus = np.zeros([60, 20])

for x in range(60):
    for y in range(20):
        matrix_modulus[x, y] = tot_mod[x][y]

fig, ax = plt.subplots(figsize=(4, 3))

ax.set_xlabel('Atom type')
ax.set_ylabel('Phonon mode')

im = ax.imshow(matrix_modulus, cmap='plasma', aspect='auto')
cbar = fig.colorbar(im, ax=ax, orientation='vertical')
cbar.set_label('Modulus phonon mode') 

ax.axvline(3.5, color='black', linestyle='--', linewidth=2)
ax.axvline(7.5, color='black', linestyle='--', linewidth=2)

ax.axhline(11.5, color='black', linewidth=2.5)

ax.set_xticks([1.5, 5.5, 13], ['S', 'Br', 'Ag'])

yticks_positions = [6, 35]
yticks_labels = ['above gap', 'below gap']

ax.set_yticks(yticks_positions)
ax.set_yticklabels(yticks_labels, rotation=90, va='center')

ax.tick_params(axis='x', which='both', length=0)
ax.tick_params(axis='y', which='both', length=0)

plt.tight_layout()
plt.savefig('modulus_phonon.pdf')