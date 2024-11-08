# Pol Benítez Colominas, November 2024
# University of Cambridge and Universitat Politècnica de Catalunya

# Plot the valence band and conduction band in the full Brillouin zone for 2D materials

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
    print(f'k-point: {kpoint[x]}, valence value: {band[x]}')

# triangle 1
cart_points_x = []
cart_points_y = []
for k in kpoint:
    cart = kpoint_to_cartesian(k, 2 * np.arctan(0.28867513 / 0.5))
    cart_points_x.append(cart[0])
    cart_points_y.append(cart[1])

triangle1 = np.array([[0, 0], kpoint_to_cartesian([0.5, 0], 2 * np.arctan(0.28867513 / 0.5)), kpoint_to_cartesian([1/3, 1/3], 2 * np.arctan(0.28867513 / 0.5))])
path1 = Path(triangle1)

grid_size = 200
xi, yi = np.linspace(-0.6, 0, grid_size), np.linspace(0, 0.6, grid_size)
xi, yi = np.meshgrid(xi, yi)

zi = griddata((cart_points_x, cart_points_y), band, (xi, yi), method='linear')

zi = np.ma.masked_invalid(zi)
mask1 = ~path1.contains_points(np.c_[xi.ravel(), yi.ravel()])

# Ensure vmin and vmax are set based on valid data points only
valid_values = zi.compressed()  # Get non-masked values
if valid_values.size > 0:
    vmin, vmax = valid_values.min(), valid_values.max()
else:
    vmin, vmax = 0, 1  # Default fallback if no valid data

heatmap = plt.pcolormesh(xi, yi, zi, shading='auto', cmap='viridis', vmin=vmin, vmax=vmax, rasterized=True)

# triangle 2
cart_points_x = []
cart_points_y = []
for k in kpoint:
    cart = kpoint_to_cartesian(k, 2 * np.arctan(0.28867513 / 0.5))
    cart_points_x.append(-cart[0])
    cart_points_y.append(cart[1])

point_triangle = kpoint_to_cartesian([1/3, 1/3], 2 * np.arctan(0.28867513 / 0.5))
point_triangle = [-point_triangle[0], point_triangle[1]]
triangle2 = np.array([[0, 0], kpoint_to_cartesian([0.5, 0], 2 * np.arctan(0.28867513 / 0.5)), point_triangle])
path2 = Path(triangle2)

grid_size = 200
xi, yi = np.linspace(0, 0.6, grid_size), np.linspace(0, 0.6, grid_size)
xi, yi = np.meshgrid(xi, yi)

zi = griddata((cart_points_x, cart_points_y), band, (xi, yi), method='linear')

zi = np.ma.masked_invalid(zi)
mask2 = ~path2.contains_points(np.c_[xi.ravel(), yi.ravel()])

heatmap = plt.pcolormesh(xi, yi, zi, shading='auto', cmap='viridis', vmin=vmin, vmax=vmax, rasterized=True)


n = 1 # rotation
# triangle 3
cart_points_x = []
cart_points_y = []
for k in kpoint:
    cart = kpoint_to_cartesian(k, 2 * np.arctan(0.28867513 / 0.5))
    cart_points_x.append(apply_rotation(cart, n*1.0472)[0])
    cart_points_y.append(apply_rotation(cart, n*1.0472)[1])

triangle3 = np.array([[0, 0], apply_rotation(kpoint_to_cartesian([0.5, 0], 2 * np.arctan(0.28867513 / 0.5)), n*1.0472), apply_rotation(kpoint_to_cartesian([1/3, 1/3], 2 * np.arctan(0.28867513 / 0.5)), n*1.0472)])
path3 = Path(triangle3)

grid_size = 200
xi, yi = np.linspace(-0.8, 0, grid_size), np.linspace(0, 0.6, grid_size)
xi, yi = np.meshgrid(xi, yi)

zi = griddata((cart_points_x, cart_points_y), band, (xi, yi), method='linear')

zi = np.ma.masked_invalid(zi)
mask3 = ~path3.contains_points(np.c_[xi.ravel(), yi.ravel()])

# Ensure vmin and vmax are set based on valid data points only
valid_values = zi.compressed()  # Get non-masked values
if valid_values.size > 0:
    vmin, vmax = valid_values.min(), valid_values.max()
else:
    vmin, vmax = 0, 1  # Default fallback if no valid data

heatmap = plt.pcolormesh(xi, yi, zi, shading='auto', cmap='viridis', vmin=vmin, vmax=vmax, rasterized=True)

# triangle 4
cart_points_x = []
cart_points_y = []
for k in kpoint:
    cart = kpoint_to_cartesian(k, 2 * np.arctan(0.28867513 / 0.5))
    cart_points_x.append(apply_rotation([-cart[0], cart[1]], n*1.0472)[0])
    cart_points_y.append(apply_rotation([-cart[0], cart[1]], n*1.0472)[1])

point_triangle = kpoint_to_cartesian([1/3, 1/3], 2 * np.arctan(0.28867513 / 0.5))
point_triangle = [-point_triangle[0], point_triangle[1]]
triangle4 = np.array([[0, 0], apply_rotation(kpoint_to_cartesian([0.5, 0], 2 * np.arctan(0.28867513 / 0.5)), n*1.0472), apply_rotation(point_triangle, n*1.0472)])
path4 = Path(triangle4)

grid_size = 200
xi, yi = np.linspace(-0.8, 0, grid_size), np.linspace(0, 0.6, grid_size)
xi, yi = np.meshgrid(xi, yi)

zi = griddata((cart_points_x, cart_points_y), band, (xi, yi), method='linear')

zi = np.ma.masked_invalid(zi)
mask4 = ~path4.contains_points(np.c_[xi.ravel(), yi.ravel()])

heatmap = plt.pcolormesh(xi, yi, zi, shading='auto', cmap='viridis', vmin=vmin, vmax=vmax, rasterized=True)


n = 2 # rotation
# triangle 5
cart_points_x = []
cart_points_y = []
for k in kpoint:
    cart = kpoint_to_cartesian(k, 2 * np.arctan(0.28867513 / 0.5))
    cart_points_x.append(apply_rotation(cart, n*1.0472)[0])
    cart_points_y.append(apply_rotation(cart, n*1.0472)[1])

triangle5 = np.array([[0, 0], apply_rotation(kpoint_to_cartesian([0.5, 0], 2 * np.arctan(0.28867513 / 0.5)), n*1.0472), apply_rotation(kpoint_to_cartesian([1/3, 1/3], 2 * np.arctan(0.28867513 / 0.5)), n*1.0472)])
path5 = Path(triangle5)

grid_size = 200
xi, yi = np.linspace(-0.8, 0, grid_size), np.linspace(-0.6, 0, grid_size)
xi, yi = np.meshgrid(xi, yi)

zi = griddata((cart_points_x, cart_points_y), band, (xi, yi), method='linear')

zi = np.ma.masked_invalid(zi)
mask5 = ~path5.contains_points(np.c_[xi.ravel(), yi.ravel()])

# Ensure vmin and vmax are set based on valid data points only
valid_values = zi.compressed()  # Get non-masked values
if valid_values.size > 0:
    vmin, vmax = valid_values.min(), valid_values.max()
else:
    vmin, vmax = 0, 1  # Default fallback if no valid data

heatmap = plt.pcolormesh(xi, yi, zi, shading='auto', cmap='viridis', vmin=vmin, vmax=vmax, rasterized=True)

# triangle 6
cart_points_x = []
cart_points_y = []
for k in kpoint:
    cart = kpoint_to_cartesian(k, 2 * np.arctan(0.28867513 / 0.5))
    cart_points_x.append(apply_rotation([-cart[0], cart[1]], n*1.0472)[0])
    cart_points_y.append(apply_rotation([-cart[0], cart[1]], n*1.0472)[1])

point_triangle = kpoint_to_cartesian([1/3, 1/3], 2 * np.arctan(0.28867513 / 0.5))
point_triangle = [-point_triangle[0], point_triangle[1]]
triangle6 = np.array([[0, 0], apply_rotation(kpoint_to_cartesian([0.5, 0], 2 * np.arctan(0.28867513 / 0.5)), n*1.0472), apply_rotation(point_triangle, n*1.0472)])
path6 = Path(triangle6)

grid_size = 200
xi, yi = np.linspace(-0.8, 0, grid_size), np.linspace(-0.6, 0, grid_size)
xi, yi = np.meshgrid(xi, yi)

zi = griddata((cart_points_x, cart_points_y), band, (xi, yi), method='linear')

zi = np.ma.masked_invalid(zi)
mask6 = ~path6.contains_points(np.c_[xi.ravel(), yi.ravel()])

heatmap = plt.pcolormesh(xi, yi, zi, shading='auto', cmap='viridis', vmin=vmin, vmax=vmax, rasterized=True)


n = 3 # rotation
# triangle 7
cart_points_x = []
cart_points_y = []
for k in kpoint:
    cart = kpoint_to_cartesian(k, 2 * np.arctan(0.28867513 / 0.5))
    cart_points_x.append(apply_rotation(cart, n*1.0472)[0])
    cart_points_y.append(apply_rotation(cart, n*1.0472)[1])

triangle7 = np.array([[0, 0], apply_rotation(kpoint_to_cartesian([0.5, 0], 2 * np.arctan(0.28867513 / 0.5)), n*1.0472), apply_rotation(kpoint_to_cartesian([1/3, 1/3], 2 * np.arctan(0.28867513 / 0.5)), n*1.0472)])
path7 = Path(triangle7)

grid_size = 200
xi, yi = np.linspace(0, 0.8, grid_size), np.linspace(-0.6, 0, grid_size)
xi, yi = np.meshgrid(xi, yi)

zi = griddata((cart_points_x, cart_points_y), band, (xi, yi), method='linear')

zi = np.ma.masked_invalid(zi)
mask7 = ~path7.contains_points(np.c_[xi.ravel(), yi.ravel()])

# Ensure vmin and vmax are set based on valid data points only
valid_values = zi.compressed()  # Get non-masked values
if valid_values.size > 0:
    vmin, vmax = valid_values.min(), valid_values.max()
else:
    vmin, vmax = 0, 1  # Default fallback if no valid data

heatmap = plt.pcolormesh(xi, yi, zi, shading='auto', cmap='viridis', vmin=vmin, vmax=vmax, rasterized=True)

# triangle 8
cart_points_x = []
cart_points_y = []
for k in kpoint:
    cart = kpoint_to_cartesian(k, 2 * np.arctan(0.28867513 / 0.5))
    cart_points_x.append(apply_rotation([-cart[0], cart[1]], n*1.0472)[0])
    cart_points_y.append(apply_rotation([-cart[0], cart[1]], n*1.0472)[1])

point_triangle = kpoint_to_cartesian([1/3, 1/3], 2 * np.arctan(0.28867513 / 0.5))
point_triangle = [-point_triangle[0], point_triangle[1]]
triangle8 = np.array([[0, 0], apply_rotation(kpoint_to_cartesian([0.5, 0], 2 * np.arctan(0.28867513 / 0.5)), n*1.0472), apply_rotation(point_triangle, n*1.0472)])
path8 = Path(triangle8)

grid_size = 200
xi, yi = np.linspace(-0.8, 0, grid_size), np.linspace(-0.6, 0, grid_size)
xi, yi = np.meshgrid(xi, yi)

zi = griddata((cart_points_x, cart_points_y), band, (xi, yi), method='linear')

zi = np.ma.masked_invalid(zi)
mask8 = ~path8.contains_points(np.c_[xi.ravel(), yi.ravel()])

heatmap = plt.pcolormesh(xi, yi, zi, shading='auto', cmap='viridis', vmin=vmin, vmax=vmax, rasterized=True)


n = 4 # rotation
# triangle 9
cart_points_x = []
cart_points_y = []
for k in kpoint:
    cart = kpoint_to_cartesian(k, 2 * np.arctan(0.28867513 / 0.5))
    cart_points_x.append(apply_rotation(cart, n*1.0472)[0])
    cart_points_y.append(apply_rotation(cart, n*1.0472)[1])

triangle9 = np.array([[0, 0], apply_rotation(kpoint_to_cartesian([0.5, 0], 2 * np.arctan(0.28867513 / 0.5)), n*1.0472), apply_rotation(kpoint_to_cartesian([1/3, 1/3], 2 * np.arctan(0.28867513 / 0.5)), n*1.0472)])
path9 = Path(triangle9)

grid_size = 200
xi, yi = np.linspace(0, 0.8, grid_size), np.linspace(-0.6, 0, grid_size)
xi, yi = np.meshgrid(xi, yi)

zi = griddata((cart_points_x, cart_points_y), band, (xi, yi), method='linear')

zi = np.ma.masked_invalid(zi)
mask9 = ~path9.contains_points(np.c_[xi.ravel(), yi.ravel()])

# Ensure vmin and vmax are set based on valid data points only
valid_values = zi.compressed()  # Get non-masked values
if valid_values.size > 0:
    vmin, vmax = valid_values.min(), valid_values.max()
else:
    vmin, vmax = 0, 1  # Default fallback if no valid data

heatmap = plt.pcolormesh(xi, yi, zi, shading='auto', cmap='viridis', vmin=vmin, vmax=vmax, rasterized=True)

# triangle 10
cart_points_x = []
cart_points_y = []
for k in kpoint:
    cart = kpoint_to_cartesian(k, 2 * np.arctan(0.28867513 / 0.5))
    cart_points_x.append(apply_rotation([-cart[0], cart[1]], n*1.0472)[0])
    cart_points_y.append(apply_rotation([-cart[0], cart[1]], n*1.0472)[1])

point_triangle = kpoint_to_cartesian([1/3, 1/3], 2 * np.arctan(0.28867513 / 0.5))
point_triangle = [-point_triangle[0], point_triangle[1]]
triangle10 = np.array([[0, 0], apply_rotation(kpoint_to_cartesian([0.5, 0], 2 * np.arctan(0.28867513 / 0.5)), n*1.0472), apply_rotation(point_triangle, n*1.0472)])
path10 = Path(triangle10)

grid_size = 200
xi, yi = np.linspace(0, 0.8, grid_size), np.linspace(-0.6, 0, grid_size)
xi, yi = np.meshgrid(xi, yi)

zi = griddata((cart_points_x, cart_points_y), band, (xi, yi), method='linear')

zi = np.ma.masked_invalid(zi)
mask10 = ~path10.contains_points(np.c_[xi.ravel(), yi.ravel()])

heatmap = plt.pcolormesh(xi, yi, zi, shading='auto', cmap='viridis', vmin=vmin, vmax=vmax, rasterized=True)


n = 5 # rotation
# triangle 11
cart_points_x = []
cart_points_y = []
for k in kpoint:
    cart = kpoint_to_cartesian(k, 2 * np.arctan(0.28867513 / 0.5))
    cart_points_x.append(apply_rotation(cart, n*1.0472)[0])
    cart_points_y.append(apply_rotation(cart, n*1.0472)[1])

triangle11 = np.array([[0, 0], apply_rotation(kpoint_to_cartesian([0.5, 0], 2 * np.arctan(0.28867513 / 0.5)), n*1.0472), apply_rotation(kpoint_to_cartesian([1/3, 1/3], 2 * np.arctan(0.28867513 / 0.5)), n*1.0472)])
path11 = Path(triangle11)

grid_size = 200
xi, yi = np.linspace(0, 0.8, grid_size), np.linspace(0, 0.6, grid_size)
xi, yi = np.meshgrid(xi, yi)

zi = griddata((cart_points_x, cart_points_y), band, (xi, yi), method='linear')

zi = np.ma.masked_invalid(zi)
mask11 = ~path11.contains_points(np.c_[xi.ravel(), yi.ravel()])

# Ensure vmin and vmax are set based on valid data points only
valid_values = zi.compressed()  # Get non-masked values
if valid_values.size > 0:
    vmin, vmax = valid_values.min(), valid_values.max()
else:
    vmin, vmax = 0, 1  # Default fallback if no valid data

heatmap = plt.pcolormesh(xi, yi, zi, shading='auto', cmap='viridis', vmin=vmin, vmax=vmax, rasterized=True)

# triangle 12
cart_points_x = []
cart_points_y = []
for k in kpoint:
    cart = kpoint_to_cartesian(k, 2 * np.arctan(0.28867513 / 0.5))
    cart_points_x.append(apply_rotation([-cart[0], cart[1]], n*1.0472)[0])
    cart_points_y.append(apply_rotation([-cart[0], cart[1]], n*1.0472)[1])

point_triangle = kpoint_to_cartesian([1/3, 1/3], 2 * np.arctan(0.28867513 / 0.5))
point_triangle = [-point_triangle[0], point_triangle[1]]
triangle12 = np.array([[0, 0], apply_rotation(kpoint_to_cartesian([0.5, 0], 2 * np.arctan(0.28867513 / 0.5)), n*1.0472), apply_rotation(point_triangle, n*1.0472)])
path12 = Path(triangle12)

grid_size = 200
xi, yi = np.linspace(0, 0.8, grid_size), np.linspace(0, 0.6, grid_size)
xi, yi = np.meshgrid(xi, yi)

zi = griddata((cart_points_x, cart_points_y), band, (xi, yi), method='linear')

zi = np.ma.masked_invalid(zi)
mask12 = ~path12.contains_points(np.c_[xi.ravel(), yi.ravel()])

heatmap = plt.pcolormesh(xi, yi, zi, shading='auto', cmap='viridis', vmin=vmin, vmax=vmax, rasterized=True)

# remove the mask
combined_mask = mask1 & mask2 & mask3 & mask4 & mask5 & mask6 & mask7 & mask8 & mask9 & mask10 & mask11 & mask12
zi.mask |= combined_mask.reshape(xi.shape)
# Color the masked regions (outside the triangle) with white
plt.imshow(combined_mask.reshape(xi.shape), extent=[-.75, .75, -0.75, 0.75], origin='lower', 
           cmap=ListedColormap(['white', 'none']))

cbar = plt.colorbar(heatmap, label='Energy (eV)') 

plt.title('Valence band TiS$_2$ monolayer')
plt.savefig('valence-band.pdf')
#################################################


############## CONDUCTION BAND ##############
plt.figure(figsize=(4, 4))

_, _, band, kpoint = get_bands('results-bands-TiS2', 9)

# shift the bands with MVB = 0
for x in range(len(band)):
    band[x] = band[x] + 4.1046
    print(f'k-point: {kpoint[x]}, conduction value: {band[x]}')

# triangle 1
cart_points_x = []
cart_points_y = []
for k in kpoint:
    cart = kpoint_to_cartesian(k, 2 * np.arctan(0.28867513 / 0.5))
    cart_points_x.append(cart[0])
    cart_points_y.append(cart[1])

triangle1 = np.array([[0, 0], kpoint_to_cartesian([0.5, 0], 2 * np.arctan(0.28867513 / 0.5)), kpoint_to_cartesian([1/3, 1/3], 2 * np.arctan(0.28867513 / 0.5))])
path1 = Path(triangle1)

grid_size = 200
xi, yi = np.linspace(-0.6, 0, grid_size), np.linspace(0, 0.6, grid_size)
xi, yi = np.meshgrid(xi, yi)

zi = griddata((cart_points_x, cart_points_y), band, (xi, yi), method='linear')

zi = np.ma.masked_invalid(zi)
mask1 = ~path1.contains_points(np.c_[xi.ravel(), yi.ravel()])

# Ensure vmin and vmax are set based on valid data points only
valid_values = zi.compressed()  # Get non-masked values
if valid_values.size > 0:
    vmin, vmax = valid_values.min(), valid_values.max()
else:
    vmin, vmax = 0, 1  # Default fallback if no valid data

heatmap = plt.pcolormesh(xi, yi, zi, shading='auto', cmap='viridis', vmin=vmin, vmax=vmax, rasterized=True)

# triangle 2
cart_points_x = []
cart_points_y = []
for k in kpoint:
    cart = kpoint_to_cartesian(k, 2 * np.arctan(0.28867513 / 0.5))
    cart_points_x.append(-cart[0])
    cart_points_y.append(cart[1])

point_triangle = kpoint_to_cartesian([1/3, 1/3], 2 * np.arctan(0.28867513 / 0.5))
point_triangle = [-point_triangle[0], point_triangle[1]]
triangle2 = np.array([[0, 0], kpoint_to_cartesian([0.5, 0], 2 * np.arctan(0.28867513 / 0.5)), point_triangle])
path2 = Path(triangle2)

grid_size = 200
xi, yi = np.linspace(0, 0.6, grid_size), np.linspace(0, 0.6, grid_size)
xi, yi = np.meshgrid(xi, yi)

zi = griddata((cart_points_x, cart_points_y), band, (xi, yi), method='linear')

zi = np.ma.masked_invalid(zi)
mask2 = ~path2.contains_points(np.c_[xi.ravel(), yi.ravel()])

heatmap = plt.pcolormesh(xi, yi, zi, shading='auto', cmap='viridis', vmin=vmin, vmax=vmax, rasterized=True)


n = 1 # rotation
# triangle 3
cart_points_x = []
cart_points_y = []
for k in kpoint:
    cart = kpoint_to_cartesian(k, 2 * np.arctan(0.28867513 / 0.5))
    cart_points_x.append(apply_rotation(cart, n*1.0472)[0])
    cart_points_y.append(apply_rotation(cart, n*1.0472)[1])

triangle3 = np.array([[0, 0], apply_rotation(kpoint_to_cartesian([0.5, 0], 2 * np.arctan(0.28867513 / 0.5)), n*1.0472), apply_rotation(kpoint_to_cartesian([1/3, 1/3], 2 * np.arctan(0.28867513 / 0.5)), n*1.0472)])
path3 = Path(triangle3)

grid_size = 200
xi, yi = np.linspace(-0.8, 0, grid_size), np.linspace(0, 0.6, grid_size)
xi, yi = np.meshgrid(xi, yi)

zi = griddata((cart_points_x, cart_points_y), band, (xi, yi), method='linear')

zi = np.ma.masked_invalid(zi)
mask3 = ~path3.contains_points(np.c_[xi.ravel(), yi.ravel()])

# Ensure vmin and vmax are set based on valid data points only
valid_values = zi.compressed()  # Get non-masked values
if valid_values.size > 0:
    vmin, vmax = valid_values.min(), valid_values.max()
else:
    vmin, vmax = 0, 1  # Default fallback if no valid data

heatmap = plt.pcolormesh(xi, yi, zi, shading='auto', cmap='viridis', vmin=vmin, vmax=vmax, rasterized=True)

# triangle 4
cart_points_x = []
cart_points_y = []
for k in kpoint:
    cart = kpoint_to_cartesian(k, 2 * np.arctan(0.28867513 / 0.5))
    cart_points_x.append(apply_rotation([-cart[0], cart[1]], n*1.0472)[0])
    cart_points_y.append(apply_rotation([-cart[0], cart[1]], n*1.0472)[1])

point_triangle = kpoint_to_cartesian([1/3, 1/3], 2 * np.arctan(0.28867513 / 0.5))
point_triangle = [-point_triangle[0], point_triangle[1]]
triangle4 = np.array([[0, 0], apply_rotation(kpoint_to_cartesian([0.5, 0], 2 * np.arctan(0.28867513 / 0.5)), n*1.0472), apply_rotation(point_triangle, n*1.0472)])
path4 = Path(triangle4)

grid_size = 200
xi, yi = np.linspace(-0.8, 0, grid_size), np.linspace(0, 0.6, grid_size)
xi, yi = np.meshgrid(xi, yi)

zi = griddata((cart_points_x, cart_points_y), band, (xi, yi), method='linear')

zi = np.ma.masked_invalid(zi)
mask4 = ~path4.contains_points(np.c_[xi.ravel(), yi.ravel()])

heatmap = plt.pcolormesh(xi, yi, zi, shading='auto', cmap='viridis', vmin=vmin, vmax=vmax, rasterized=True)


n = 2 # rotation
# triangle 5
cart_points_x = []
cart_points_y = []
for k in kpoint:
    cart = kpoint_to_cartesian(k, 2 * np.arctan(0.28867513 / 0.5))
    cart_points_x.append(apply_rotation(cart, n*1.0472)[0])
    cart_points_y.append(apply_rotation(cart, n*1.0472)[1])

triangle5 = np.array([[0, 0], apply_rotation(kpoint_to_cartesian([0.5, 0], 2 * np.arctan(0.28867513 / 0.5)), n*1.0472), apply_rotation(kpoint_to_cartesian([1/3, 1/3], 2 * np.arctan(0.28867513 / 0.5)), n*1.0472)])
path5 = Path(triangle5)

grid_size = 200
xi, yi = np.linspace(-0.8, 0, grid_size), np.linspace(-0.6, 0, grid_size)
xi, yi = np.meshgrid(xi, yi)

zi = griddata((cart_points_x, cart_points_y), band, (xi, yi), method='linear')

zi = np.ma.masked_invalid(zi)
mask5 = ~path5.contains_points(np.c_[xi.ravel(), yi.ravel()])

# Ensure vmin and vmax are set based on valid data points only
valid_values = zi.compressed()  # Get non-masked values
if valid_values.size > 0:
    vmin, vmax = valid_values.min(), valid_values.max()
else:
    vmin, vmax = 0, 1  # Default fallback if no valid data

heatmap = plt.pcolormesh(xi, yi, zi, shading='auto', cmap='viridis', vmin=vmin, vmax=vmax, rasterized=True)

# triangle 6
cart_points_x = []
cart_points_y = []
for k in kpoint:
    cart = kpoint_to_cartesian(k, 2 * np.arctan(0.28867513 / 0.5))
    cart_points_x.append(apply_rotation([-cart[0], cart[1]], n*1.0472)[0])
    cart_points_y.append(apply_rotation([-cart[0], cart[1]], n*1.0472)[1])

point_triangle = kpoint_to_cartesian([1/3, 1/3], 2 * np.arctan(0.28867513 / 0.5))
point_triangle = [-point_triangle[0], point_triangle[1]]
triangle6 = np.array([[0, 0], apply_rotation(kpoint_to_cartesian([0.5, 0], 2 * np.arctan(0.28867513 / 0.5)), n*1.0472), apply_rotation(point_triangle, n*1.0472)])
path6 = Path(triangle6)

grid_size = 200
xi, yi = np.linspace(-0.8, 0, grid_size), np.linspace(-0.6, 0, grid_size)
xi, yi = np.meshgrid(xi, yi)

zi = griddata((cart_points_x, cart_points_y), band, (xi, yi), method='linear')

zi = np.ma.masked_invalid(zi)
mask6 = ~path6.contains_points(np.c_[xi.ravel(), yi.ravel()])

heatmap = plt.pcolormesh(xi, yi, zi, shading='auto', cmap='viridis', vmin=vmin, vmax=vmax, rasterized=True)


n = 3 # rotation
# triangle 7
cart_points_x = []
cart_points_y = []
for k in kpoint:
    cart = kpoint_to_cartesian(k, 2 * np.arctan(0.28867513 / 0.5))
    cart_points_x.append(apply_rotation(cart, n*1.0472)[0])
    cart_points_y.append(apply_rotation(cart, n*1.0472)[1])

triangle7 = np.array([[0, 0], apply_rotation(kpoint_to_cartesian([0.5, 0], 2 * np.arctan(0.28867513 / 0.5)), n*1.0472), apply_rotation(kpoint_to_cartesian([1/3, 1/3], 2 * np.arctan(0.28867513 / 0.5)), n*1.0472)])
path7 = Path(triangle7)

grid_size = 200
xi, yi = np.linspace(0, 0.8, grid_size), np.linspace(-0.6, 0, grid_size)
xi, yi = np.meshgrid(xi, yi)

zi = griddata((cart_points_x, cart_points_y), band, (xi, yi), method='linear')

zi = np.ma.masked_invalid(zi)
mask7 = ~path7.contains_points(np.c_[xi.ravel(), yi.ravel()])

# Ensure vmin and vmax are set based on valid data points only
valid_values = zi.compressed()  # Get non-masked values
if valid_values.size > 0:
    vmin, vmax = valid_values.min(), valid_values.max()
else:
    vmin, vmax = 0, 1  # Default fallback if no valid data

heatmap = plt.pcolormesh(xi, yi, zi, shading='auto', cmap='viridis', vmin=vmin, vmax=vmax, rasterized=True)

# triangle 8
cart_points_x = []
cart_points_y = []
for k in kpoint:
    cart = kpoint_to_cartesian(k, 2 * np.arctan(0.28867513 / 0.5))
    cart_points_x.append(apply_rotation([-cart[0], cart[1]], n*1.0472)[0])
    cart_points_y.append(apply_rotation([-cart[0], cart[1]], n*1.0472)[1])

point_triangle = kpoint_to_cartesian([1/3, 1/3], 2 * np.arctan(0.28867513 / 0.5))
point_triangle = [-point_triangle[0], point_triangle[1]]
triangle8 = np.array([[0, 0], apply_rotation(kpoint_to_cartesian([0.5, 0], 2 * np.arctan(0.28867513 / 0.5)), n*1.0472), apply_rotation(point_triangle, n*1.0472)])
path8 = Path(triangle8)

grid_size = 200
xi, yi = np.linspace(-0.8, 0, grid_size), np.linspace(-0.6, 0, grid_size)
xi, yi = np.meshgrid(xi, yi)

zi = griddata((cart_points_x, cart_points_y), band, (xi, yi), method='linear')

zi = np.ma.masked_invalid(zi)
mask8 = ~path8.contains_points(np.c_[xi.ravel(), yi.ravel()])

heatmap = plt.pcolormesh(xi, yi, zi, shading='auto', cmap='viridis', vmin=vmin, vmax=vmax, rasterized=True)


n = 4 # rotation
# triangle 9
cart_points_x = []
cart_points_y = []
for k in kpoint:
    cart = kpoint_to_cartesian(k, 2 * np.arctan(0.28867513 / 0.5))
    cart_points_x.append(apply_rotation(cart, n*1.0472)[0])
    cart_points_y.append(apply_rotation(cart, n*1.0472)[1])

triangle9 = np.array([[0, 0], apply_rotation(kpoint_to_cartesian([0.5, 0], 2 * np.arctan(0.28867513 / 0.5)), n*1.0472), apply_rotation(kpoint_to_cartesian([1/3, 1/3], 2 * np.arctan(0.28867513 / 0.5)), n*1.0472)])
path9 = Path(triangle9)

grid_size = 200
xi, yi = np.linspace(0, 0.8, grid_size), np.linspace(-0.6, 0, grid_size)
xi, yi = np.meshgrid(xi, yi)

zi = griddata((cart_points_x, cart_points_y), band, (xi, yi), method='linear')

zi = np.ma.masked_invalid(zi)
mask9 = ~path9.contains_points(np.c_[xi.ravel(), yi.ravel()])

# Ensure vmin and vmax are set based on valid data points only
valid_values = zi.compressed()  # Get non-masked values
if valid_values.size > 0:
    vmin, vmax = valid_values.min(), valid_values.max()
else:
    vmin, vmax = 0, 1  # Default fallback if no valid data

heatmap = plt.pcolormesh(xi, yi, zi, shading='auto', cmap='viridis', vmin=vmin, vmax=vmax, rasterized=True)

# triangle 10
cart_points_x = []
cart_points_y = []
for k in kpoint:
    cart = kpoint_to_cartesian(k, 2 * np.arctan(0.28867513 / 0.5))
    cart_points_x.append(apply_rotation([-cart[0], cart[1]], n*1.0472)[0])
    cart_points_y.append(apply_rotation([-cart[0], cart[1]], n*1.0472)[1])

point_triangle = kpoint_to_cartesian([1/3, 1/3], 2 * np.arctan(0.28867513 / 0.5))
point_triangle = [-point_triangle[0], point_triangle[1]]
triangle10 = np.array([[0, 0], apply_rotation(kpoint_to_cartesian([0.5, 0], 2 * np.arctan(0.28867513 / 0.5)), n*1.0472), apply_rotation(point_triangle, n*1.0472)])
path10 = Path(triangle10)

grid_size = 200
xi, yi = np.linspace(0, 0.8, grid_size), np.linspace(-0.6, 0, grid_size)
xi, yi = np.meshgrid(xi, yi)

zi = griddata((cart_points_x, cart_points_y), band, (xi, yi), method='linear')

zi = np.ma.masked_invalid(zi)
mask10 = ~path10.contains_points(np.c_[xi.ravel(), yi.ravel()])

heatmap = plt.pcolormesh(xi, yi, zi, shading='auto', cmap='viridis', vmin=vmin, vmax=vmax, rasterized=True)


n = 5 # rotation
# triangle 11
cart_points_x = []
cart_points_y = []
for k in kpoint:
    cart = kpoint_to_cartesian(k, 2 * np.arctan(0.28867513 / 0.5))
    cart_points_x.append(apply_rotation(cart, n*1.0472)[0])
    cart_points_y.append(apply_rotation(cart, n*1.0472)[1])

triangle11 = np.array([[0, 0], apply_rotation(kpoint_to_cartesian([0.5, 0], 2 * np.arctan(0.28867513 / 0.5)), n*1.0472), apply_rotation(kpoint_to_cartesian([1/3, 1/3], 2 * np.arctan(0.28867513 / 0.5)), n*1.0472)])
path11 = Path(triangle11)

grid_size = 200
xi, yi = np.linspace(0, 0.8, grid_size), np.linspace(0, 0.6, grid_size)
xi, yi = np.meshgrid(xi, yi)

zi = griddata((cart_points_x, cart_points_y), band, (xi, yi), method='linear')

zi = np.ma.masked_invalid(zi)
mask11 = ~path11.contains_points(np.c_[xi.ravel(), yi.ravel()])

# Ensure vmin and vmax are set based on valid data points only
valid_values = zi.compressed()  # Get non-masked values
if valid_values.size > 0:
    vmin, vmax = valid_values.min(), valid_values.max()
else:
    vmin, vmax = 0, 1  # Default fallback if no valid data

heatmap = plt.pcolormesh(xi, yi, zi, shading='auto', cmap='viridis', vmin=vmin, vmax=vmax, rasterized=True)

# triangle 12
cart_points_x = []
cart_points_y = []
for k in kpoint:
    cart = kpoint_to_cartesian(k, 2 * np.arctan(0.28867513 / 0.5))
    cart_points_x.append(apply_rotation([-cart[0], cart[1]], n*1.0472)[0])
    cart_points_y.append(apply_rotation([-cart[0], cart[1]], n*1.0472)[1])

point_triangle = kpoint_to_cartesian([1/3, 1/3], 2 * np.arctan(0.28867513 / 0.5))
point_triangle = [-point_triangle[0], point_triangle[1]]
triangle12 = np.array([[0, 0], apply_rotation(kpoint_to_cartesian([0.5, 0], 2 * np.arctan(0.28867513 / 0.5)), n*1.0472), apply_rotation(point_triangle, n*1.0472)])
path12 = Path(triangle12)

grid_size = 200
xi, yi = np.linspace(0, 0.8, grid_size), np.linspace(0, 0.6, grid_size)
xi, yi = np.meshgrid(xi, yi)

zi = griddata((cart_points_x, cart_points_y), band, (xi, yi), method='linear')

zi = np.ma.masked_invalid(zi)
mask12 = ~path12.contains_points(np.c_[xi.ravel(), yi.ravel()])

heatmap = plt.pcolormesh(xi, yi, zi, shading='auto', cmap='viridis', vmin=vmin, vmax=vmax, rasterized=True)

# remove the mask
combined_mask = mask1 & mask2 & mask3 & mask4 & mask5 & mask6 & mask7 & mask8 & mask9 & mask10 & mask11 & mask12
zi.mask |= combined_mask.reshape(xi.shape)
# Color the masked regions (outside the triangle) with white
plt.imshow(combined_mask.reshape(xi.shape), extent=[-.75, .75, -0.75, 0.75], origin='lower', 
           cmap=ListedColormap(['white', 'none']))

cbar = plt.colorbar(heatmap, label='Energy (eV)') 

plt.title('Conduction band TiS$_2$ monolayer')
plt.savefig('conduction-band.pdf')
#################################################
