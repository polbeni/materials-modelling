import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.ticker import MultipleLocator

##### total #####
def compute_color(df, band, index):
    """
    It computes the RGB color plus alpha of a given band and k-point and the energy

    Inputs:
        df -> dataframe containing the data
        band -> band we want to plot
        index -> index of the k-point
    """

    k_point = df.loc[index, 'k_length']

    energy = df.loc[index, 'band'+ str(band)]

    Ag = df.loc[index, 'Wan3InBand' + str(band)] + df.loc[index, 'Wan4InBand' + str(band)] + df.loc[index, 'Wan5InBand' + str(band)] + df.loc[index, 'Wan6InBand' + str(band)] + df.loc[index, 'Wan7InBand' + str(band)] + df.loc[index, 'Wan8InBand' + str(band)] + df.loc[index, 'Wan9InBand' + str(band)] + df.loc[index, 'Wan10InBand' + str(band)] + df.loc[index, 'Wan11InBand' + str(band)] + df.loc[index, 'Wan12InBand' + str(band)] + df.loc[index, 'Wan13InBand' + str(band)] + df.loc[index, 'Wan14InBand' + str(band)] + df.loc[index, 'Wan15InBand' + str(band)] + df.loc[index, 'Wan16InBand' + str(band)] + df.loc[index, 'Wan17InBand' + str(band)] + df.loc[index, 'Wan18InBand' + str(band)] + df.loc[index, 'Wan19InBand' + str(band)] + df.loc[index, 'Wan20InBand' + str(band)] + df.loc[index, 'Wan31InBand' + str(band)] + df.loc[index, 'Wan32InBand' + str(band)] + df.loc[index, 'Wan33InBand' + str(band)] + df.loc[index, 'Wan34InBand' + str(band)] + df.loc[index, 'Wan35InBand' + str(band)] + df.loc[index, 'Wan36InBand' + str(band)] + df.loc[index, 'Wan37InBand' + str(band)] + df.loc[index, 'Wan38InBand' + str(band)] + df.loc[index, 'Wan39InBand' + str(band)] 

    S = df.loc[index, 'Wan0InBand' + str(band)] + df.loc[index, 'Wan1InBand' + str(band)] + df.loc[index, 'Wan2InBand' + str(band)] + df.loc[index, 'Wan25InBand' + str(band)] + df.loc[index, 'Wan26InBand' + str(band)] + df.loc[index, 'Wan27InBand' + str(band)] + df.loc[index, 'Wan28InBand' + str(band)] + df.loc[index, 'Wan29InBand' + str(band)] + df.loc[index, 'Wan30InBand' + str(band)]

    Br = df.loc[index, 'Wan21InBand' + str(band)] + df.loc[index, 'Wan22InBand' + str(band)] + df.loc[index, 'Wan23InBand' + str(band)] + df.loc[index, 'Wan24InBand' + str(band)] + df.loc[index, 'Wan40InBand' + str(band)] + df.loc[index, 'Wan41InBand' + str(band)] + df.loc[index, 'Wan42InBand' + str(band)] + df.loc[index, 'Wan43InBand' + str(band)] + df.loc[index, 'Wan44InBand' + str(band)]
   
    R = Br / (Ag + S + Br)
    G = Ag / (Ag + S + Br)
    B = S / (Ag + S + Br)

    return k_point, energy, R, G, B





fig, ax = plt.subplots(1, 2, figsize=(8, 3))

ax[0].set_title('Non-distorted')

ax[0].set_xlim(0, 2.278)
ax[0].set_ylim(-2, 4)
ax[0].set_ylabel('Energy (eV)')

ax[0].axvline(0.34992, color='black', linestyle='--', linewidth=.8)
ax[0].axvline(0.693648, color='black', linestyle='--', linewidth=.8)
ax[0].axvline(1.18413, color='black', linestyle='--', linewidth=.8)
ax[0].axvline(1.785, color='black', linestyle='--', linewidth=.8)

ax[0].set_xticks(ticks=[0, 0.34992, 0.693648, 1.18413, 1.785, 2.278], labels=['$\\Gamma$', 'X', 'M', '$\\Gamma$', 'R', 'X'])



df = pd.read_csv('data/BSAg3SBr_orig.txt', delim_whitespace=True)


for band_num in range(14, 28):
    print(band_num)
    k_val = []
    energy = []
    color = []
    for x in range(1004):
        k, e, r, g, b = compute_color(df,band_num, x)
        k_val.append(k)
        energy.append(e)
        color.append([r, g, b])

    points = np.array([k_val, energy]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    lc = LineCollection(segments, colors=color, linewidth=2)

    ax[0].add_collection(lc)

ax[0].plot([100, 102], [100, 100], color=(0, 1, 0), label='Ag')
ax[0].plot([100, 102], [100, 100], color=(0, 0, 1), label='S')
ax[0].plot([100, 102], [100, 100], color=(1, 0, 0), label='Br')
ax[0].legend(facecolor='white', framealpha=1)

major_locator = MultipleLocator(1)
minor_locator = MultipleLocator(.2)
ax[0].yaxis.set_major_locator(major_locator)
ax[0].yaxis.set_minor_locator(minor_locator)



ax[1].set_title('Distorted')

ax[1].set_xlim(0, 2.278)
ax[1].set_ylim(-2, 4)
#ax[1].set_ylabel('$E-E_{F}$ (eV)')

ax[1].axvline(0.34992, color='black', linestyle='--', linewidth=.8)
ax[1].axvline(0.693648, color='black', linestyle='--', linewidth=.8)
ax[1].axvline(1.18413, color='black', linestyle='--', linewidth=.8)
ax[1].axvline(1.785, color='black', linestyle='--', linewidth=.8)

ax[1].set_xticks(ticks=[0, 0.34992, 0.693648, 1.18413, 1.785, 2.278], labels=['$\\Gamma$', 'X', 'M', '$\\Gamma$', 'R', 'X'])





df = pd.read_csv('data/BSAg3SBr_dis.txt', delim_whitespace=True)


for band_num in range(14, 28):
    print(band_num)
    k_val = []
    energy = []
    color = []
    for x in range(1004):
        k, e, r, g, b = compute_color(df,band_num, x)
        k_val.append(k)
        energy.append(e)
        color.append([r, g, b])

    points = np.array([k_val, energy]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    lc = LineCollection(segments, colors=color, linewidth=2)

    ax[1].add_collection(lc)







ax[1].plot([100, 102], [100, 100], color=(0, 1, 0), label='Ag')
ax[1].plot([100, 102], [100, 100], color=(0, 0, 1), label='S')
ax[1].plot([100, 102], [100, 100], color=(1, 0, 0), label='Br')
ax[1].legend(facecolor='white', framealpha=1)

major_locator = MultipleLocator(1)
minor_locator = MultipleLocator(.2)
ax[1].yaxis.set_major_locator(major_locator)
ax[1].yaxis.set_minor_locator(minor_locator)


plt.tight_layout()
plt.savefig('bands-tight-binding.pdf')
