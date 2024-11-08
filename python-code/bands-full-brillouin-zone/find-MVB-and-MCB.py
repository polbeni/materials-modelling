# Pol Benítez Colominas, November 2024
# University of Cambridge and Universitat Politècnica de Catalunya

# Plot near maximum and minimum of valence band and conduction band in the symmetric zone of the Brillouin zone (for a 2D material)

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.colors import ListedColormap
from scipy.interpolate import griddata

def get_bands(OUTCAR_path, nband):
    """
    Get the energy bands from the OUTCAR file from a VASP calculation for a given band

    Inputs:
        OUTCAR_path -> path to the OUTCAR file
        nband - > number of the band we want the results (in the order they appear in the OUTCAR file)
    Outputs:
        fermi_energy -> value of the fermi energy
        num_kpoints -> number of k-points in the calculation
        band_values -> array with the value of the band for each k-point
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

    band_values = np.zeros(num_kpoints)
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

            if (band + 1) == nband:
                band_values[kpoint] = float(line.split()[1])
        
        line = OUTCAR.readline()

    OUTCAR.close()

    return fermi_energy, num_kpoints, band_values, k_points

def kpoint_to_cartesian(k_point, angle):
    """
    It express a k-point in cartesian coordinates and rotates the desired angle

    Inputs:
        k_point -> coordinates of the k-point
        angle -> angle we want to rotate the points
    """

    x_coord = (k_point[0] * 1) + (k_point[1] * 0)
    y_coord = (k_point[0] * (1 / np.sqrt(3))) + (k_point[1] * (2 / np.sqrt(3)))


    final_x_coord = (x_coord * np.cos(angle)) - (y_coord * np.sin(angle))
    final_y_coord = (x_coord * np.sin(angle)) + (y_coord * np.cos(angle))

    cart_kpoint = np.array([final_x_coord, final_y_coord])

    return cart_kpoint

def apply_rotation(point, angle):
    """
    It express a k-point in cartesian coordinates and rotates the desired angle

    Inputs:
        point -> cartesian coordinates of the point
        angle -> angle we want to rotate the points
    """


    final_x_coord = (point[0] * np.cos(angle)) - (point[1] * np.sin(angle))
    final_y_coord = (point[0] * np.sin(angle)) + (point[1] * np.cos(angle))

    final = np.array([final_x_coord, final_y_coord])

    return final


############## VALENCE BAND ##############
plt.figure(figsize=(4, 4))

_, _, band, kpoint = get_bands('results-bands-TiS2', 8)

# shift the bands with MVB = 0
for x in range(len(band)):
    band[x] = band[x] + 4.1046

n = 9

cart_points_x = []
cart_points_y = []
for k in kpoint:
    cart = kpoint_to_cartesian(k, 2 * np.arctan(0.28867513 / 0.5) + n*1.0472)
    cart_points_x.append(cart[0])
    cart_points_y.append(cart[1])

triangle = np.array([[0, 0], kpoint_to_cartesian([0.5, 0], 2 * np.arctan(0.28867513 / 0.5) + n*1.0472), kpoint_to_cartesian([1/3, 1/3], 2 * np.arctan(0.28867513 / 0.5) + n*1.0472)])
path = Path(triangle)

grid_size = 200
xi, yi = np.linspace(0, 0.6, grid_size), np.linspace(-0.6, 0, grid_size)
xi, yi = np.meshgrid(xi, yi)

zi = griddata((cart_points_x, cart_points_y), band, (xi, yi), method='linear')

zi = np.ma.masked_invalid(zi)
zi = np.ma.masked_invalid(zi)
mask = ~path.contains_points(np.c_[xi.ravel(), yi.ravel()])
zi.mask |= mask.reshape(xi.shape)

# Define the z-value range you want to display
z_min, z_max = -0.1, 0  # Adjust these values as needed

# Create an additional mask based on the z-value range
range_mask = (zi < z_min) | (zi > z_max)
zi.mask |= range_mask

# Ensure vmin and vmax are set based on valid data points only
valid_values = zi.compressed()  # Get non-masked values
if valid_values.size > 0:
    vmin, vmax = valid_values.min(), valid_values.max()
else:
    vmin, vmax = 0, 1  # Default fallback if no valid data

heatmap = plt.pcolormesh(xi, yi, zi, shading='auto', cmap='viridis', vmin=vmin, vmax=vmax, rasterized=True)


# Color the masked regions (outside the triangle) with white
plt.imshow(mask.reshape(xi.shape), extent=[-0.0666, .4, -0.65, 0.08], origin='lower', 
           cmap=ListedColormap(['white', 'none']))

cbar = plt.colorbar(heatmap, label='Energy (eV)') 

path1_x = [0, kpoint_to_cartesian([0.5, 0], 2 * np.arctan(0.28867513 / 0.5) + n*1.0472)[0]]
path1_y = [0, kpoint_to_cartesian([0.5, 0], 2 * np.arctan(0.28867513 / 0.5) + n*1.0472)[1]]

path2_x = [0, kpoint_to_cartesian([1/3, 1/3], 2 * np.arctan(0.28867513 / 0.5) + n*1.0472)[0]]
path2_y = [0, kpoint_to_cartesian([1/3, 1/3], 2 * np.arctan(0.28867513 / 0.5) + n*1.0472)[1]]

path3_x = [kpoint_to_cartesian([0.5, 0], 2 * np.arctan(0.28867513 / 0.5) + n*1.0472)[0], kpoint_to_cartesian([1/3, 1/3], 2 * np.arctan(0.28867513 / 0.5) + n*1.0472)[0]]
path3_y = [kpoint_to_cartesian([0.5, 0], 2 * np.arctan(0.28867513 / 0.5) + n*1.0472)[1], kpoint_to_cartesian([1/3, 1/3], 2 * np.arctan(0.28867513 / 0.5) + n*1.0472)[1]]

plt.plot(path1_x, path1_y, color='black', linewidth=1.5)
plt.plot(path2_x, path2_y, color='black', linewidth=1.5)
plt.plot(path3_x, path3_y, color='black', linewidth=1.5)

plt.text(-0.02, 0.02, '$\\Gamma$', fontsize=15)
plt.text(-0.02, -0.63, 'M', fontsize=15)
plt.text(0.32, -0.63, 'K', fontsize=15)
plt.text(0.19, -0.31, '$\\Lambda$', fontsize=15)
plt.text(-0.04, -0.31, '$\\Sigma$', fontsize=15)

plt.xlim(-0.05, 0.4)

plt.title('Valence band TiS$_2$ monolayer')
plt.savefig('valence-band-zoom-max.pdf')
plt.savefig('valence-band-zoom-max.png', dpi=500)
#################################################


############## CONDUCTION BAND ##############
plt.figure(figsize=(4, 4))

_, _, band, kpoint = get_bands('results-bands-TiS2', 9)

# shift the bands with MVB = 0
for x in range(len(band)):
    band[x] = band[x] + 4.1046

n = 9

cart_points_x = []
cart_points_y = []
for k in kpoint:
    cart = kpoint_to_cartesian(k, 2 * np.arctan(0.28867513 / 0.5) + n*1.0472)
    cart_points_x.append(cart[0])
    cart_points_y.append(cart[1])

triangle = np.array([[0, 0], kpoint_to_cartesian([0.5, 0], 2 * np.arctan(0.28867513 / 0.5) + n*1.0472), kpoint_to_cartesian([1/3, 1/3], 2 * np.arctan(0.28867513 / 0.5) + n*1.0472)])
path = Path(triangle)

grid_size = 200
xi, yi = np.linspace(0, 0.6, grid_size), np.linspace(-0.6, 0, grid_size)
xi, yi = np.meshgrid(xi, yi)

zi = griddata((cart_points_x, cart_points_y), band, (xi, yi), method='linear')

zi = np.ma.masked_invalid(zi)
zi = np.ma.masked_invalid(zi)
mask = ~path.contains_points(np.c_[xi.ravel(), yi.ravel()])
zi.mask |= mask.reshape(xi.shape)

# Define the z-value range you want to display
z_min, z_max = 1.5, 1.7  # Adjust these values as needed

# Create an additional mask based on the z-value range
range_mask = (zi < z_min) | (zi > z_max)
zi.mask |= range_mask

# Ensure vmin and vmax are set based on valid data points only
valid_values = zi.compressed()  # Get non-masked values
if valid_values.size > 0:
    vmin, vmax = valid_values.min(), valid_values.max()
else:
    vmin, vmax = 0, 1  # Default fallback if no valid data

heatmap = plt.pcolormesh(xi, yi, zi, shading='auto', cmap='viridis', vmin=vmin, vmax=vmax, rasterized=True)


# Color the masked regions (outside the triangle) with white
plt.imshow(mask.reshape(xi.shape), extent=[-0.0666, .4, -0.65, 0.08], origin='lower', 
           cmap=ListedColormap(['white', 'none']))

cbar = plt.colorbar(heatmap, label='Energy (eV)') 

path1_x = [0, kpoint_to_cartesian([0.5, 0], 2 * np.arctan(0.28867513 / 0.5) + n*1.0472)[0]]
path1_y = [0, kpoint_to_cartesian([0.5, 0], 2 * np.arctan(0.28867513 / 0.5) + n*1.0472)[1]]

path2_x = [0, kpoint_to_cartesian([1/3, 1/3], 2 * np.arctan(0.28867513 / 0.5) + n*1.0472)[0]]
path2_y = [0, kpoint_to_cartesian([1/3, 1/3], 2 * np.arctan(0.28867513 / 0.5) + n*1.0472)[1]]

path3_x = [kpoint_to_cartesian([0.5, 0], 2 * np.arctan(0.28867513 / 0.5) + n*1.0472)[0], kpoint_to_cartesian([1/3, 1/3], 2 * np.arctan(0.28867513 / 0.5) + n*1.0472)[0]]
path3_y = [kpoint_to_cartesian([0.5, 0], 2 * np.arctan(0.28867513 / 0.5) + n*1.0472)[1], kpoint_to_cartesian([1/3, 1/3], 2 * np.arctan(0.28867513 / 0.5) + n*1.0472)[1]]

plt.plot(path1_x, path1_y, color='black', linewidth=1.5)
plt.plot(path2_x, path2_y, color='black', linewidth=1.5)
plt.plot(path3_x, path3_y, color='black', linewidth=1.5)

plt.text(-0.02, 0.02, '$\\Gamma$', fontsize=15)
plt.text(-0.02, -0.63, 'M', fontsize=15)
plt.text(0.32, -0.63, 'K', fontsize=15)
plt.text(0.19, -0.31, '$\\Lambda$', fontsize=15)
plt.text(-0.04, -0.31, '$\\Sigma$', fontsize=15)

plt.xlim(-0.05, 0.4)

plt.title('Conduction band TiS$_2$ monolayer')
plt.savefig('conduction-band-zoom-min.pdf')
plt.savefig('conduction-band-zoom-min.png', dpi=500)
#################################################