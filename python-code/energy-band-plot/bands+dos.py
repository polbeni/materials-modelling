# Pol Benítez Colominas, September-October 2024
# University of Cambridge and Universitat Politècnica de Catalunya

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.gridspec import GridSpec
from scipy.interpolate import interp1d


####################################################### IMPORTANT #######################################################
# give the path to the OUTCAR file
path_results = 'data/OUTCAR' # path to OUTCAR

# define here the k-path looking at the OUTCAR or KPOITNS files, the number in the k_points_in_path correspond 
# to their position in the files, starting from 0. Then give a name to these points.
k_points_in_path = [[0, 1, 2, 3], [3, 6, 8, 9], [9, 7, 4, 0], [0, 10, 16, 19], [19, 17, 12, 3]] # index of k-points of interest
sym_points_in_path = [0, 3, 9, 0, 19, 3] # index of high symmetry points in the path
names_sym_points = ['$\\Gamma$', 'X', 'M', '$\\Gamma$', 'R', 'X'] # names of the high symmetry points

# define the density of your k_point mesh
mesh_density = 1000 # number of points in the interpolated k-path

# set the energy limit values to represent in the plot
min_energy = -7.5 # min value of the energy band to represent in the plot
max_energy = 12.5 # max value of the energy band to represent in the plot
#########################################################################################################################


# functions
def get_bands(OUTCAR_path):
    """
    Get the energy bands from the OUTCAR file from a VASP calculation

    Inputs:
        OUTCAR_path -> path to the OUTCAR file
    Outputs:
        fermi_energy -> value of the fermi energy
        num_kpoints -> number of k-points in the calculation
        num_bands -> total number of energy bands
        bands -> matrix with the energy bands with dimension num_kpoints x num_bands
        k_points -> an array with the position of each k-point in the calculation
    """

    with open(OUTCAR_path, 'r') as file:
        searc_pattern = 'E-fermi'

        for line in file:
            if searc_pattern in line:
                fermi_energy = float(line.split()[2])

    with open(OUTCAR_path, 'r') as file:
        searc_pattern = 'NKPTS'

        for line in file:
            if searc_pattern in line:
                num_kpoints = int(line.split()[3])
                num_bands = int(line.split()[14])

    total_num_lines = 0
    finish_count = False
    with open(OUTCAR_path, 'r') as file:
        searc_pattern = 'E-fermi'

        for line in file:
            if finish_count == False:
                total_num_lines = total_num_lines + 1
            if searc_pattern in line:
                finish_count = True

    bands = np.zeros([num_kpoints, num_bands])
    k_points = []

    OUTCAR = open(OUTCAR_path, 'r')

    for x in range(total_num_lines + 3):
        line = OUTCAR.readline()

    for kpoint in range(num_kpoints):
        line = OUTCAR.readline()

        k_pos = np.array([float(line.split()[3]), float(line.split()[4]), float(line.split()[5])])
        k_points.append(k_pos)

        line = OUTCAR.readline()

        for band in range(num_bands):
            line = OUTCAR.readline()

            bands[kpoint, band] = float(line.split()[1])
        
        line = OUTCAR.readline()

    OUTCAR.close()

    return fermi_energy, num_kpoints, num_bands, bands, k_points

def bands_segment(electronic_bands, segments, n_band):
    """
    It gives the values of the bands for each segment in the k-path

    Inputs:
        electronic_bands -> values of the band for all the path
        segments -> array with different arrays representing each segment in the k-path
        n_bands -> number of the electronic band
    Outputs:
        bands -> array with arrays of the energy bands for the different segments in the k-path
        band_points -> an array with all the energy in k-path (without duplicate values)
    """

    bands = []

    for segment in segments:
        band = []
        for k_point in segment:
            band.append(electronic_bands[k_point, n_band])

        bands.append(band)

    it_seg = 0
    band_points = []
    for segment in bands:
        if it_seg == 0:
            for elem in segment:
                band_points.append(elem)
        else:
            it_el = 0
            for elem in segment:
                if it_el != 0:
                    band_points.append(elem)
                it_el = it_el + 1
        it_seg = it_seg + 1

    return bands, band_points

def find_minCB_maxVB(fermi_energy, bands):
    """
    Finds the minimum of the valence band and the maximum of the conduction band (for semiconductors and insulators)

    Inputs:
        fermi_energy -> energy value of the fermi level
        bands -> values for all the energy bands
    Outputs:
        minCB -> minimum of the conduction band
        maxVB -> maximum of the valence band
        bandgap -> band gap of the material
        type_of_gap -> if the material is a direct or indirect semiconductor
    """

    minCB = 1e5
    maxVB = -1e5

    k_VB = 0
    k_CB = 0

    for band in bands:
        it_el = 0
        for element in band:
            if (abs(element - fermi_energy) < abs(minCB - fermi_energy)) and (abs(element - fermi_energy) > 5e-2) and (element > fermi_energy):
                minCB = element
                k_CB = it_el
            if (abs(element - fermi_energy) < abs(maxVB - fermi_energy)) and (abs(element - fermi_energy) > 5e-2) and (element < fermi_energy):
                maxVB = element
                k_VB = it_el

            it_el = it_el + 1

    bandgap = minCB - maxVB

    if k_CB == k_VB:
        type_of_gap = 'Direct'
    else:
        type_of_gap = 'Indirect'

    return minCB, maxVB, bandgap, type_of_gap


# get the values
fermi_energy, n_kpoints, n_bands, bands, k_points = get_bands(path_results)

# interpolate the bands
distance_between_sym_points = []
for point in range(len(sym_points_in_path) - 1):
    distance = np.linalg.norm(k_points[sym_points_in_path[point + 1]] - k_points[sym_points_in_path[point]])
    distance_between_sym_points.append(distance)

points_to_interpol = []
total_num_k_interp = 0
for element in distance_between_sym_points:
    points_to_interpol.append(int(element * mesh_density))
    
    total_num_k_interp = total_num_k_interp + int(element * mesh_density)

final_k_path = np.linspace(0, 1, total_num_k_interp)

k_new_interval_norm = []
num_points = 0
initial = 0
for segment in points_to_interpol:
    num_points = num_points + segment
    final = num_points / total_num_k_interp

    k_new_interval_norm.append(np.linspace(initial, final, segment))

    initial = final

total_distance = sum(distance_between_sym_points)
k_old_interval_norm = []
initial = 0
distance = 0
it_loop = 0
for element in distance_between_sym_points:
    distance = distance + element
    final = distance / total_distance

    k_old_interval_norm.append(np.linspace(initial, final, len(k_points_in_path[it_loop])))

    initial = final
    it_loop = it_loop + 1

total_num_k_interp = total_num_k_interp - len(points_to_interpol) + 1
k_points_paths_interp = []

band_matrix_interp = []

total_old_segment = []
it_segment = 0
for segment in k_old_interval_norm:
    if it_segment == 0:
        for element in segment:
            total_old_segment.append(element)
    else:
        it_element = 0
        for element in segment:
            if it_element != 0:
                total_old_segment.append(element)

            it_element = it_element + 1
    it_segment = it_segment + 1

final_k_interval = np.linspace(0, 1, mesh_density)
for band in range(n_bands):
    _, band_segment = bands_segment(bands, k_points_in_path, band)

    interpolated_band = interp1d(total_old_segment, band_segment, kind='quadratic', fill_value='extrapolate')(final_k_interval)

    band_matrix_interp.append(interpolated_band)
        

minCB, maxVB, bg, type_gap = find_minCB_maxVB(fermi_energy, band_matrix_interp)
print('################### ELECTRONIC RESULTS ###################')
print(f'           The value of the band gap is:  {bg:.3f} eV')
print(f'        The gap of the semiconductor is:  {type_gap}')
print(f'              The fermi level energy is: {fermi_energy:.3f} eV')
print(f'  The minimum of the conduction band is: {minCB:.3f} eV')
print(f'     The maximum of the valence band is: {maxVB:.3f} eV')
print('##########################################################')

# adjust the bands to fermi_energy=0
band_matrix_interp_shift = []
for band in band_matrix_interp:
    new_band = []
    for element in band:
        new_band.append(element - fermi_energy)
    band_matrix_interp_shift.append(new_band)
        

# plot the bands for the desired path
plt.figure()
fig, axs = plt.subplots(1, 2, figsize=(4, 2))

gs = GridSpec(1, 2, width_ratios=[1, 0.3])

axs[0] = plt.subplot(gs[0,0])
axs[1] = plt.subplot(gs[0,1])

it_segment = 0
for segment in k_old_interval_norm:
    if it_segment != 0:
        axs[0].axvline(segment[0], color='black', linestyle='--', linewidth=0.8)
    
    it_segment = it_segment + 1

axs[0].set_xlim(0, 1)
axs[0].set_ylim(min_energy, max_energy)

axs[0].set_ylabel('$E-E_{F}$ (eV)')

for band in range(n_bands):
    if band < 23:
        axs[0].plot(final_k_interval, band_matrix_interp_shift[band], color='firebrick', alpha=0.8)
    else:
        axs[0].plot(final_k_interval, band_matrix_interp_shift[band], color='royalblue', alpha=0.8)

new_labels = names_sym_points
old_labels = []
it_segment = 0
for segment in k_new_interval_norm:
    if it_segment == 0:
        old_labels.append(segment[0])

    old_labels.append(segment[-1])

    it_segment = it_segment + 1
axs[0].set_xticks(ticks=old_labels, labels=new_labels)


axs[0].plot(0.5204268292682926, minCB - fermi_energy, 'o', markersize=8, color='navy')
axs[0].plot(0.3048780487804878, maxVB - fermi_energy, 'o', markersize=8, color='darkred')

major_locator = MultipleLocator(5)
minor_locator = MultipleLocator(1)
axs[0].yaxis.set_major_locator(major_locator)
axs[0].yaxis.set_minor_locator(minor_locator)


# eDOS
path = 'data/TDOS.dat'
total_dos = np.loadtxt(path)
path = 'data/PDOS_Ag.dat'
Ag_dos = np.loadtxt(path)
path = 'data/PDOS_S.dat'
S_dos = np.loadtxt(path)
path = 'data/PDOS_Br.dat'
Br_dos = np.loadtxt(path)

axs[1].plot(total_dos[:,-1], total_dos[:,0], color='black', linewidth=.9, label='Total')
axs[1].fill_between(total_dos[:,-1], total_dos[:,0], color='grey', edgecolor='None', alpha=0.5)
axs[1].plot(Ag_dos[:,-1], Ag_dos[:,0], color='cornflowerblue', linewidth=.9, label='Ag')
axs[1].fill_between(Ag_dos[:,-1], Ag_dos[:,0], color='cornflowerblue', edgecolor='None', alpha=0.5)
axs[1].plot(S_dos[:,-1], S_dos[:,0], color='lime', linewidth=.9, label='S')
axs[1].fill_between(S_dos[:,-1], S_dos[:,0], color='lime', edgecolor='None', alpha=0.5)
axs[1].plot(Br_dos[:,-1], Br_dos[:,0], color='violet', linewidth=.9, label='Br')
axs[1].fill_between(Br_dos[:,-1], Br_dos[:,0], color='violet', edgecolor='None', alpha=0.5)
axs[1].set_xlim(0,6)
axs[1].set_ylim(-7.5, 12.5)
axs[1].set_yticklabels([])

major_locator = MultipleLocator(5)
minor_locator = MultipleLocator(1)
axs[1].yaxis.set_major_locator(major_locator)
axs[1].yaxis.set_minor_locator(minor_locator)

major_locator = MultipleLocator(3)
minor_locator = MultipleLocator(0.6)
axs[1].xaxis.set_major_locator(major_locator)
axs[1].xaxis.set_minor_locator(minor_locator)

#axs[1].legend(frameon=False, ncol=2, fontsize=6)

plt.tight_layout()
plt.savefig('bands.pdf')