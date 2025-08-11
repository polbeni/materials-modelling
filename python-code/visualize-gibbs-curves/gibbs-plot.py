# Pol Benítez Colominas, August 2025
# Universitat Politècnica de Catalunya

# Visualize the Gibbs curves with respect (p,T)

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib.colors as mcolors
from scipy.interpolate import interp1d
from mpl_toolkits.mplot3d import Axes3D

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


def check_valid_T(pressure_in, temperature_in, material):
    """
    It verifies the minimum T at what the given p computations are valid (no extrapolated)
    It looks for all the phases at takes the maximum minimum T of the phases

    Inputs: 
        pressure_in -> desired pressure (its index)
        temperature_in -> desired temperature (its index)
    """
    
    # define condition, True if there is a volume smaller, False otherwise
    condition = False 

    if material == 'monoGS':
        val_monoGS = read_volume(f'data-pressure/monoGS/volume-temperature-{pressure_in + 1}.dat', temperature_in)
        if val_monoGS <= 132.80:
            condition = True

    if material == 'cubic':
        val_cubic = read_volume(f'data-pressure/cubic/volume-temperature-{pressure_in + 1}.dat', temperature_in)
        if val_cubic <= 123.90:
            condition = True

    if material == 'orthoAPref':
        val_orthoAPref = read_volume(f'data-pressure/orthoAPref/volume-temperature-{pressure_in + 1}.dat', temperature_in)
        if val_orthoAPref <= 141.68:
            condition = True

    if material == 'orthoI':
        val_orthoI = read_volume(f'data-pressure/orthoI/volume-temperature-{pressure_in + 1}.dat', temperature_in)
        if val_orthoI <= 254.50:
            condition = True

    if material == 'orthoIIIFE':
        val_orthoIIIFE = read_volume(f'data-pressure/orthoIIIFE/volume-temperature-{pressure_in + 1}.dat', temperature_in)
        if val_orthoIIIFE <= 128.00:
            condition = True

    if material == 'orthoIstar':
        val_orthoIstar = read_volume(f'data-pressure/orthoIstar/volume-temperature-{pressure_in + 1}.dat', temperature_in)
        if val_orthoIstar <= 266.20:
            condition = True

    if material == 'tetra':
        val_tetra = read_volume(f'data-pressure/tetra/volume-temperature-{pressure_in + 1}.dat', temperature_in)
        if val_tetra <= 126.70:
            condition = True

    return condition


############ ANHARMONIC ############

pressure_points = 200
temperature_points = 181

pressure_array = np.linspace(0, 10, pressure_points) # in GPa
temperature_array = np.linspace(0, 1800, temperature_points) # in K

energy_orthoIstar = np.zeros([temperature_points, pressure_points])
energy_orthoI = np.zeros([temperature_points, pressure_points])
energy_orthoIIIFE = np.zeros([temperature_points, pressure_points])
energy_orthoAPref = np.zeros([temperature_points, pressure_points])
energy_tetra = np.zeros([temperature_points, pressure_points])
energy_cubic = np.zeros([temperature_points, pressure_points])

for pressure in range(pressure_points):
    for temperature in range(temperature_points):
        condition = check_valid_T(pressure, temperature, 'monoGS')
        if condition == False:
            _, gibbs = get_gibbs_energy(f'data-pressure/monoGS/gibbs-{pressure + 1}.dat', 12)
            energy_monoGS = gibbs[temperature] * 1e3

            condition = check_valid_T(pressure, temperature, 'orthoIstar')
            if condition == False:
                _, gibbs = get_gibbs_energy(f'data-pressure/orthoIstar/gibbs-{pressure + 1}.dat', 24)
                energy_orthoIstar[temperature, pressure] = (gibbs[temperature]*1e3 - energy_monoGS)
            else:
                energy_orthoIstar[temperature, pressure] = np.nan

            condition = check_valid_T(pressure, temperature, 'orthoI')
            if condition == False:
                _, gibbs = get_gibbs_energy(f'data-pressure/orthoI/gibbs-{pressure + 1}.dat', 24)
                energy_orthoI[temperature, pressure] = (gibbs[temperature]*1e3 - energy_monoGS)
            else:
                energy_orthoI[temperature, pressure] = np.nan

            condition = check_valid_T(pressure, temperature, 'orthoIIIFE')
            if condition == False:
                _, gibbs = get_gibbs_energy(f'data-pressure/orthoIIIFE/gibbs-{pressure + 1}.dat', 12)
                energy_orthoIIIFE[temperature, pressure] = (gibbs[temperature]*1e3 - energy_monoGS)
            else:
                energy_orthoIIIFE[temperature, pressure] = np.nan

            condition = check_valid_T(pressure, temperature, 'orthoAPref')
            if condition == False:
                _, gibbs = get_gibbs_energy(f'data-pressure/orthoAPref/gibbs-{pressure + 1}.dat', 12)
                energy_orthoAPref[temperature, pressure] = (gibbs[temperature]*1e3 - energy_monoGS)
            else:
                energy_orthoAPref[temperature, pressure] = np.nan

            condition = check_valid_T(pressure, temperature, 'tetra')
            if condition == False:
                _, gibbs = get_gibbs_energy(f'data-pressure/tetra/gibbs-{pressure + 1}.dat', 12)
                energy_tetra[temperature, pressure] = (gibbs[temperature]*1e3 - energy_monoGS)
            else:
                energy_tetra[temperature, pressure] = np.nan

            condition = check_valid_T(pressure, temperature, 'cubic')
            if condition == False:
                _, gibbs = get_gibbs_energy(f'data-pressure/cubic/gibbs-{pressure + 1}.dat', 12)
                energy_cubic[temperature, pressure] = (gibbs[temperature]*1e3 - energy_monoGS)
            else:
                energy_cubic[temperature, pressure] = np.nan
        else:
            energy_orthoIstar[temperature, pressure] = np.nan
            energy_orthoI[temperature, pressure] = np.nan
            energy_orthoIIIFE[temperature, pressure] = np.nan
            energy_orthoAPref[temperature, pressure] = np.nan
            energy_tetra[temperature, pressure] = np.nan
            energy_cubic[temperature, pressure] = np.nan

P, T = np.meshgrid(pressure_array, temperature_array)

fig = plt.figure(figsize=(4, 4))
ax = fig.add_subplot(111, projection='3d')
ax.set_box_aspect([1, 1, 1.4])

# Create a z=0 plane
z0 = np.zeros_like(P)
ax.plot_surface(P, T, z0, color='gray', alpha=0.7, edgecolor='none')

alpha_gibbs_surfaces = 0.6

ax.plot_surface(P, T, energy_orthoI, color='lightsteelblue', alpha=alpha_gibbs_surfaces)
ax.plot_surface(P, T, energy_orthoIstar, color='lightcoral', alpha=alpha_gibbs_surfaces)
ax.plot_surface(P, T, energy_orthoIIIFE, color='seagreen', alpha=alpha_gibbs_surfaces)
ax.plot_surface(P, T, energy_orthoAPref, color='rebeccapurple', alpha=alpha_gibbs_surfaces)
ax.plot_surface(P, T, energy_tetra, color='burlywood', alpha=alpha_gibbs_surfaces)
ax.plot_surface(P, T, energy_cubic, color='cornflowerblue', alpha=alpha_gibbs_surfaces)

ax.set_xlim(0,10)
ax.set_ylim(0, 1800)
ax.set_zlim(-70, 215)

ax.set_xlabel('Pressure (GPa)')
ax.set_ylabel('Temperature (K)')
zlabel = ax.set_zlabel('$G_{P2_1/c}-G_{x}$ (meV)')

zlabel.set_rotation(90)  # Try 90 or 270 for vertical alignment
ax.zaxis.set_rotate_label(False)  # Disable auto-rotation
zlabel.set_verticalalignment('bottom')  # Adjust vertical position
zlabel.set_horizontalalignment('center')

ax.view_init(elev=10, azim=-140)

ax.set_yticks(np.arange(0, 1801, 400))


# Define the corners of the box
x = [0, 10]
y = [0, 1800]
z = [-70, 215]

# List of all 12 edges
edges = [
    # bottom square
    ([x[0], x[1]], [y[0], y[0]], [z[0], z[0]]),
    ([x[0], x[1]], [y[1], y[1]], [z[0], z[0]]),
    ([x[0], x[0]], [y[0], y[1]], [z[0], z[0]]),
    ([x[1], x[1]], [y[0], y[1]], [z[0], z[0]]),

    # top square
    ([x[0], x[1]], [y[0], y[0]], [z[1], z[1]]),
    ([x[0], x[1]], [y[1], y[1]], [z[1], z[1]]),
    ([x[0], x[0]], [y[0], y[1]], [z[1], z[1]]),
    ([x[1], x[1]], [y[0], y[1]], [z[1], z[1]]),

    # vertical edges
    ([x[0], x[0]], [y[0], y[0]], [z[0], z[1]]),
    ([x[0], x[0]], [y[1], y[1]], [z[0], z[1]]),
    ([x[1], x[1]], [y[0], y[0]], [z[0], z[1]]),
    ([x[1], x[1]], [y[1], y[1]], [z[0], z[1]])
]

linewidth_box = 1

# Plot all edges
for edge in edges:
    ax.plot(*edge, color='black', linewidth=linewidth_box)

ax.plot([5, 5], [y[0], y[0]], [z[0], z[1]], color='black', linewidth=linewidth_box)
ax.plot([5, 5], [y[1], y[1]], [z[0], z[1]], color='black', linewidth=linewidth_box)
ax.plot([5, 5], [y[0], y[1]], [z[0], z[0]], color='black', linewidth=linewidth_box)
ax.plot([5, 5], [y[0], y[1]], [z[1], z[1]], color='black', linewidth=linewidth_box)


def print_intersection(surface1, surface2, color_line):
    """
    It prints the intersection between two surfaces

    Inputs:
        surface1, surface2: the two surfaces
        color_line: color of the intersection line
    """
    
    Z1 = surface1
    Z2 = surface2

    diff = Z1 - Z2

    threshold = 0.1  # in meV, adjust as needed
    intersection_mask = np.abs(diff) < threshold

    P_line = P[intersection_mask]
    T_line = T[intersection_mask]
    Z_line = Z1[intersection_mask]

    ax.plot3D(P_line, T_line, Z_line, color=color_line, linewidth=2)

# metastable to stable
print_intersection(z0, energy_orthoI, 'red')
print_intersection(energy_orthoI, energy_tetra, 'red')

# metastable to metastable
print_intersection(energy_cubic, energy_orthoI, 'blue')
print_intersection(energy_orthoAPref, energy_orthoI, 'blue')
print_intersection(energy_orthoIIIFE, energy_orthoI, 'blue')
print_intersection(energy_orthoIstar, energy_orthoI, 'blue')

print_intersection(energy_orthoIstar, energy_cubic, 'blue')
print_intersection(energy_orthoIIIFE, energy_cubic, 'blue')
print_intersection(energy_orthoAPref, energy_cubic, 'blue')

print_intersection(energy_orthoIIIFE, energy_orthoIstar, 'blue')
print_intersection(energy_orthoAPref, energy_orthoIstar, 'blue')

print_intersection(energy_orthoAPref, energy_orthoIIIFE, 'blue')

print_intersection(energy_orthoAPref, energy_tetra, 'blue')
print_intersection(energy_cubic, energy_tetra, 'blue')
print_intersection(energy_orthoIIIFE, energy_tetra, 'blue')
print_intersection(energy_orthoIstar, energy_tetra, 'blue')



plt.tight_layout()
plt.savefig('gibbs_3d.pdf', bbox_inches='tight', pad_inches=0.2)
#plt.show()


# plot the isobar surfaces
def print_intersection_curve(curve1, curve2, color_point):
    """
    It prints the intersection between two curves

    Inputs:
        surface1, surface2: the two surfaces
        color_line: color of the intersection point
    """
    
    Z1 = curve1
    Z2 = curve2

    diff = Z1 - Z2

    threshold = 0.45  # in meV, adjust as needed
    intersection_mask = []
    condition = False
    for x in range(len(diff)):
        if condition == False:
            if np.abs(diff[x]) < threshold:
                intersection_mask.append(True)
                condition = True
            else:
                intersection_mask.append(False)
        else:
            intersection_mask.append(False)

    T_line = T[intersection_mask]
    Z_line = Z1[intersection_mask]

    axs.plot(T_line[:, 0], Z_line, marker='o', linestyle='', color=color_point)


fig, axs = plt.subplots(figsize=(4, 2.5))

axs.set_title('$p=0$ GPa')
axs.set_xlabel('Temperature (K)')
axs.set_ylabel('$G_{P2_1/c}-G_{x}$ (meV)')

axs.set_xlim(0, 1800)
axs.set_ylim(-70, 215)

axs.axhline(0, color='gray', linewidth=0.7, alpha=0.7)

axs.plot(T[:, 0], energy_orthoI[:, 0], color='lightsteelblue', alpha=alpha_gibbs_surfaces)
axs.plot(T[:, 0], energy_orthoIstar[0:, 0], color='lightcoral', alpha=alpha_gibbs_surfaces)
axs.plot(T[0:160, 0], energy_orthoIIIFE[0:160, 0], color='seagreen', alpha=alpha_gibbs_surfaces)
axs.plot(T[:, 0], energy_orthoAPref[:, 0], color='rebeccapurple', alpha=alpha_gibbs_surfaces)
axs.plot(T[:, 0], energy_tetra[:, 0], color='burlywood', alpha=alpha_gibbs_surfaces)
axs.plot(T[:, 0], energy_cubic[:, 0], color='cornflowerblue', alpha=alpha_gibbs_surfaces)

major_locator = MultipleLocator(400)
minor_locator = MultipleLocator(80)
axs.xaxis.set_major_locator(major_locator)
axs.xaxis.set_minor_locator(minor_locator)

major_locator = MultipleLocator(50)
minor_locator = MultipleLocator(10)
axs.yaxis.set_major_locator(major_locator)
axs.yaxis.set_minor_locator(minor_locator)

print_intersection_curve(energy_orthoI[:, 0], energy_orthoIstar[:, 0], 'blue')
print_intersection_curve(energy_orthoIIIFE[:, 0], energy_orthoIstar[:, 0], 'blue')
print_intersection_curve(energy_orthoAPref[:, 0], energy_orthoIstar[:, 0], 'blue')
print_intersection_curve(energy_tetra[:, 0], energy_orthoIstar[:, 0], 'blue')


plt.tight_layout()
plt.savefig('isobar_0.pdf')



fig, axs = plt.subplots(figsize=(4, 2.5))

axs.set_title('$p=5$ GPa')
axs.set_xlabel('Temperature (K)')
axs.set_ylabel('$G_{P2_1/c}-G_{x}$ (meV)')

axs.set_xlim(0, 1800)
axs.set_ylim(-70, 215)

axs.axhline(0, color='gray', linewidth=0.7, alpha=0.7)

axs.plot(T[:, 0], energy_orthoI[:, int(len(pressure_array)/2)], color='lightsteelblue', alpha=alpha_gibbs_surfaces)
axs.plot(T[:, 0], energy_orthoIstar[:, int(len(pressure_array)/2)], color='lightcoral', alpha=alpha_gibbs_surfaces)
axs.plot(T[:, 0], energy_orthoIIIFE[:, int(len(pressure_array)/2)], color='seagreen', alpha=alpha_gibbs_surfaces)
axs.plot(T[:, 0], energy_orthoAPref[:, int(len(pressure_array)/2)], color='rebeccapurple', alpha=alpha_gibbs_surfaces)
axs.plot(T[:, 0], energy_tetra[:, int(len(pressure_array)/2)], color='burlywood', alpha=alpha_gibbs_surfaces)
axs.plot(T[:, 0], energy_cubic[:, int(len(pressure_array)/2)], color='cornflowerblue', alpha=alpha_gibbs_surfaces)

major_locator = MultipleLocator(400)
minor_locator = MultipleLocator(80)
axs.xaxis.set_major_locator(major_locator)
axs.xaxis.set_minor_locator(minor_locator)

major_locator = MultipleLocator(50)
minor_locator = MultipleLocator(10)
axs.yaxis.set_major_locator(major_locator)
axs.yaxis.set_minor_locator(minor_locator)

print_intersection_curve(energy_tetra[:, int(len(pressure_array)/2)], [0]*len(temperature_array), 'red')

print_intersection_curve(energy_tetra[:, int(len(pressure_array)/2)], energy_orthoIstar[:, int(len(pressure_array)/2)], 'blue')
print_intersection_curve(energy_tetra[:, int(len(pressure_array)/2)], energy_orthoIIIFE[:, int(len(pressure_array)/2)], 'blue')
print_intersection_curve(energy_orthoAPref[:, int(len(pressure_array)/2)], energy_cubic[:, int(len(pressure_array)/2)], 'blue')

plt.tight_layout()
plt.savefig('isobar_5.pdf')



fig, axs = plt.subplots(figsize=(4, 2.5))

axs.set_title('$p=10$ GPa')
axs.set_xlabel('Temperature (K)')
axs.set_ylabel('$G_{P2_1/c}-G_{x}$ (meV)')

axs.set_xlim(0, 1800)
axs.set_ylim(-70, 215)

axs.axhline(0, color='gray', linewidth=0.7, alpha=0.7)

axs.plot(T[:, 0], energy_orthoI[:, -1], color='lightsteelblue', alpha=alpha_gibbs_surfaces)
axs.plot(T[:, 0], energy_orthoIstar[:, -1], color='lightcoral', alpha=alpha_gibbs_surfaces)
axs.plot(T[:, 0], energy_orthoIIIFE[:, -1], color='seagreen', alpha=alpha_gibbs_surfaces)
axs.plot(T[:, 0], energy_orthoAPref[:, -1], color='rebeccapurple', alpha=alpha_gibbs_surfaces)
axs.plot(T[:, 0], energy_tetra[:, -1], color='burlywood', alpha=alpha_gibbs_surfaces)
axs.plot(T[:, 0], energy_cubic[:, -1], color='cornflowerblue', alpha=alpha_gibbs_surfaces)

major_locator = MultipleLocator(400)
minor_locator = MultipleLocator(80)
axs.xaxis.set_major_locator(major_locator)
axs.xaxis.set_minor_locator(minor_locator)

major_locator = MultipleLocator(50)
minor_locator = MultipleLocator(10)
axs.yaxis.set_major_locator(major_locator)
axs.yaxis.set_minor_locator(minor_locator)

plt.tight_layout()
plt.savefig('isobar_10.pdf')


# labels of materials
fig, axs = plt.subplots(figsize=(4, 2.5))

axs.set_title('$p=10$ GPa')
axs.set_xlabel('Temperature (K)')
axs.set_ylabel('$G_{P2_1/c}-G_{x}$ (meV)')

axs.set_xlim(0, 1800)
axs.set_ylim(-70, 215)

axs.axhline(0, color='gray', linewidth=0.7, alpha=0.7)

axs.plot(T[:, 0], energy_orthoI[:, 0], color='lightsteelblue', alpha=alpha_gibbs_surfaces, label='$Pbca$ I')
axs.plot(T[:, 0], energy_orthoIstar[:, 0], color='lightcoral', alpha=alpha_gibbs_surfaces, label='$Pbca$ II')
axs.plot(T[:, 0], energy_orthoIIIFE[:, 0], color='seagreen', alpha=alpha_gibbs_surfaces, label='$Pca2_1$')
axs.plot(T[:, 0], energy_orthoAPref[:, 0], color='rebeccapurple', alpha=alpha_gibbs_surfaces, label='$Pbcn$')
axs.plot(T[:, 0], energy_tetra[:, 0], color='burlywood', alpha=alpha_gibbs_surfaces, label='$P4_2/nmc$')
axs.plot(T[:, 0], energy_cubic[:, 0], color='cornflowerblue', alpha=alpha_gibbs_surfaces, label='$Fm\overline{3}m$')

major_locator = MultipleLocator(400)
minor_locator = MultipleLocator(80)
axs.xaxis.set_major_locator(major_locator)
axs.xaxis.set_minor_locator(minor_locator)

major_locator = MultipleLocator(50)
minor_locator = MultipleLocator(10)
axs.yaxis.set_major_locator(major_locator)
axs.yaxis.set_minor_locator(minor_locator)

axs.legend(bbox_to_anchor=(1, 0.8), frameon=False)
plt.tight_layout()
plt.savefig('labels.pdf')
