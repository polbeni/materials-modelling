# Pol Benítez Colominas, August 2025
# Universitat Politècnica de Catalunya

# Show the movement of ions from AIMD calculation

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib.ticker import MultipleLocator

from mpl_toolkits.mplot3d.art3d import Line3DCollection
from matplotlib import cm
from matplotlib.colors import Normalize, LinearSegmentedColormap

def get_movement(path_XDATCAR, ion_number, number_steps, nblocks, time_step):
    """
    This function returns a matrix with the times and positions for a given atom
    from a XDATCAR file generated with VASP ab intio molecular dynamics simulations

    Inputs:
        path_XDATCAR: the path of our XDATCAR file
        ion_number: the position in the POSCAR file of the atom we want to study
        number_steps: number of steps that we want to read in the XDATCAR file
        nblocks: interval of configurations that we take
        time_step: the time step of our simulation
    """

    XDATCAR_file = open(path_XDATCAR, "r")

    for _ in range(2):
        actual_line = XDATCAR_file.readline()

    lattice_parameters = []
    actual_line = XDATCAR_file.readline()
    lattice_parameters.append([float(actual_line.split()[0]), float(actual_line.split()[1]), float(actual_line.split()[2])])
    actual_line = XDATCAR_file.readline()
    lattice_parameters.append([float(actual_line.split()[0]), float(actual_line.split()[1]), float(actual_line.split()[2])])
    actual_line = XDATCAR_file.readline()
    lattice_parameters.append([float(actual_line.split()[0]), float(actual_line.split()[1]), float(actual_line.split()[2])])

    for _ in range(2):
        actual_line = XDATCAR_file.readline()
    
    num_atoms = int(actual_line.split()[0]) + int(actual_line.split()[1])

    ion_movement = np.zeros((int(number_steps/nblocks), 4))

    num_line = num_atoms 
    line_counter = num_atoms 
    nblocks_counter = nblocks
    num_configuration = 0

    num_iterations = (num_atoms + 1)*number_steps
    for iteration in range(num_iterations):
        actual_line = XDATCAR_file.readline()
        

        if actual_line.split()[0] == 'Direct':
            nblocks_counter = nblocks_counter + 1
        else:
            line_counter = line_counter + 1
            if (((line_counter - ion_number) % num_line) == 0) and ((nblocks_counter % nblocks) == 0):
                ion_movement[num_configuration, 0] = nblocks_counter*time_step
                ion_movement[num_configuration, 1] = (float(actual_line.split()[0]) * lattice_parameters[0][0]) + (float(actual_line.split()[1]) * lattice_parameters[1][0]) + (float(actual_line.split()[2]) * lattice_parameters[2][0])
                ion_movement[num_configuration, 2] = (float(actual_line.split()[0]) * lattice_parameters[0][1]) + (float(actual_line.split()[1]) * lattice_parameters[1][1]) + (float(actual_line.split()[2]) * lattice_parameters[2][1])
                ion_movement[num_configuration, 3] = (float(actual_line.split()[0]) * lattice_parameters[0][2]) + (float(actual_line.split()[1]) * lattice_parameters[1][2]) + (float(actual_line.split()[2]) * lattice_parameters[2][2])

                num_configuration = num_configuration + 1
    
    XDATCAR_file.close()

    return ion_movement


def plot_traj(path_XDATCAR, ion_number, number_steps, nblocks, time_step, color_code):
    """
    Plots the trajectroy of the desired atom

    Inputs:
        path_XDATCAR: the path of our XDATCAR file
        ion_number: the position in the POSCAR file of the atom we want to study
        number_steps: number of steps that we want to read in the XDATCAR file
        nblocks: interval of configurations that we take
        time_step: the time step of our simulation
        color_code: starting and final color for the color gradient to show the temporal change
    """

    ion_m = get_movement(path_XDATCAR, ion_number, number_steps, nblocks, time_step)

    #ion_m[:,1] = ion_m[:,1] - ion_m[0,1]
    #ion_m[:,2] = ion_m[:,2] - ion_m[0,2]
    #ion_m[:,3] = ion_m[:,3] - ion_m[0,3]

    points = ion_m[:, 1:].reshape(-1, 1, 3)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    custom_cmap = LinearSegmentedColormap.from_list(f"{color_code[0]}_{color_code[1]}", color_code)

    norm = Normalize(ion_m[:,0].min(), ion_m[:,0].max())

    colors = custom_cmap(norm(ion_m[:-1, 0]))

    lc = Line3DCollection(segments, colors=colors, linewidth=1.5, alpha=0.7)
    ax.add_collection3d(lc)

"""

####### CUBIC 300K #######
fig = plt.figure(figsize=(4, 4))
ax = fig.add_subplot(111, projection='3d')

ax.set_title('$T=300$ K')

path = 'data/cubic/300K/XDATCAR'
number_steps = 30000
nblocks = 10
time_step = 0.0015

atom_subcell = 4 # to denote the subcell in the supercell

atoms_Hf = []
atom = 1 + atom_subcell
for _ in range(4):
    atoms_Hf.append(atom)
    atom = atom + 27

atoms_O = []
atom = 109 + atom_subcell
for _ in range(8):
    atoms_O.append(atom)
    atom = atom + 27

for atom in atoms_Hf:
    plot_traj(path, atom, number_steps, nblocks, time_step, ['brown', 'green'])

for atom in atoms_O:
    plot_traj(path, atom, number_steps, nblocks, time_step, ['red', 'yellow'])


ax.set_xlabel('$x$ (Å)')
ax.set_ylabel('$y$ (Å)')
zlabel = ax.set_zlabel('$z$ (Å)')
zlabel.set_rotation(90)  # Try 90 or 270 for vertical alignment
ax.zaxis.set_rotate_label(False)  # Disable auto-rotation
zlabel.set_verticalalignment('bottom')  # Adjust vertical position
zlabel.set_horizontalalignment('center')

ax.view_init(elev=30, azim=225)

plt.tight_layout()
plt.savefig('cubic-300.pdf', bbox_inches='tight', pad_inches=0.35)
#plt.show()
plt.close()
##########################


####### CUBIC 2500K #######
fig = plt.figure(figsize=(4, 4))
ax = fig.add_subplot(111, projection='3d')

ax.set_title('$T=2500$ K')

path = 'data/cubic/2500K/XDATCAR'
number_steps = 30000
nblocks = 10
time_step = 0.0015

atom_subcell = 4 # to denote the subcell in the supercell

atoms_Hf = []
atom = 1 + atom_subcell
for _ in range(4):
    atoms_Hf.append(atom)
    atom = atom + 27

atoms_O = []
atom = 109 + atom_subcell
for _ in range(8):
    atoms_O.append(atom)
    atom = atom + 27

for atom in atoms_Hf:
    plot_traj(path, atom, number_steps, nblocks, time_step, ['brown', 'green'])

for atom in atoms_O:
    plot_traj(path, atom, number_steps, nblocks, time_step, ['red', 'yellow'])

ax.set_xlabel('$x$ (Å)')
ax.set_ylabel('$y$ (Å)')
zlabel = ax.set_zlabel('$z$ (Å)')
zlabel.set_rotation(90)  # Try 90 or 270 for vertical alignment
ax.zaxis.set_rotate_label(False)  # Disable auto-rotation
zlabel.set_verticalalignment('bottom')  # Adjust vertical position
zlabel.set_horizontalalignment('center')

ax.view_init(elev=30, azim=225)

plt.tight_layout()
plt.savefig('cubic-2500.pdf', bbox_inches='tight', pad_inches=0.35)
#plt.show()
plt.close()
##########################



##### CUBIC diffusing ####
fig = plt.figure(figsize=(4, 4))
ax = fig.add_subplot(111, projection='3d')

path = 'data/cubic/300K/XDATCAR'
number_steps = 30000
nblocks = 10
time_step = 0.0015

atom_subcell = 4 # to denote the subcell in the supercell

atoms_O = []
atom = 109 + atom_subcell
for _ in range(8):
    atoms_O.append(atom)
    atom = atom + 27

plot_traj(path, atoms_O[4], number_steps, nblocks, time_step, ['red', 'yellow'])

ax.view_init(elev=30, azim=225)

plt.tight_layout()
plt.savefig('cubic-300-O.pdf', bbox_inches='tight', pad_inches=0.35)
#plt.show()
plt.close()

fig = plt.figure(figsize=(4, 4))
ax = fig.add_subplot(111, projection='3d')

path = 'data/cubic/2500K/XDATCAR'
number_steps = 30000
nblocks = 10
time_step = 0.0015

atom_subcell = 4 # to denote the subcell in the supercell

atoms_O = []
atom = 109 + atom_subcell
for _ in range(8):
    atoms_O.append(atom)
    atom = atom + 27

plot_traj(path, atoms_O[4], number_steps, nblocks, time_step, ['red', 'yellow'])

ax.view_init(elev=30, azim=225)

plt.tight_layout()
plt.savefig('cubic-2500-O.pdf', bbox_inches='tight', pad_inches=0.35)
#plt.show()
plt.close()




fig = plt.figure(figsize=(4, 4))
ax = fig.add_subplot(111, projection='3d')

path = 'data/cubic/300K/XDATCAR'
number_steps = 30000
nblocks = 10
time_step = 0.0015

atom_subcell = 4 # to denote the subcell in the supercell

atoms_Hf = []
atom = 1 + atom_subcell
for _ in range(4):
    atoms_Hf.append(atom)
    atom = atom + 27

plot_traj(path, atoms_Hf[0], number_steps, nblocks, time_step, ['brown', 'green'])

ax.view_init(elev=30, azim=225)

plt.tight_layout()
plt.savefig('cubic-300-Hf.pdf', bbox_inches='tight', pad_inches=0.35)
#plt.show()
plt.close()

fig = plt.figure(figsize=(4, 4))
ax = fig.add_subplot(111, projection='3d')

path = 'data/cubic/2500K/XDATCAR'
number_steps = 30000
nblocks = 10
time_step = 0.0015

atom_subcell = 4 # to denote the subcell in the supercell

atoms_Hf = []
atom = 1 + atom_subcell
for _ in range(4):
    atoms_Hf.append(atom)
    atom = atom + 27

plot_traj(path, atoms_Hf[0], number_steps, nblocks, time_step, ['brown', 'green'])

ax.view_init(elev=30, azim=225)

plt.tight_layout()
plt.savefig('cubic-2500-Hf.pdf', bbox_inches='tight', pad_inches=0.35)
#plt.show()
plt.close()
##########################

"""


####### MONOIII 300K #######
fig = plt.figure(figsize=(4, 4))
ax = fig.add_subplot(111, projection='3d')

ax.set_title('$T=300$ K')

path = 'data/monoIII/300K/XDATCAR'
number_steps = 30000
nblocks = 10
time_step = 0.0015

atom_subcell = 0 # to denote the subcell in the supercell

atoms_Hf = []
atom = 1 + atom_subcell
for _ in range(8):
    atoms_Hf.append(atom)
    
    ion_m = get_movement(path, atom, number_steps, nblocks, time_step)
    print(ion_m[0,1], ion_m[0,2], ion_m[0,3])

    atom = atom + 18

atoms_O = []
atom = 145 + atom_subcell
for _ in range(16):
    atoms_O.append(atom)

    ion_m = get_movement(path, atom, number_steps, nblocks, time_step)
    print(ion_m[0,1], ion_m[0,2], ion_m[0,3])

    atom = atom + 18

for atom in atoms_Hf:
    plot_traj(path, atom, number_steps, nblocks, time_step, ['brown', 'green'])

for atom in atoms_O:
    plot_traj(path, atom, number_steps, nblocks, time_step, ['red', 'yellow'])

ax.set_xlabel('$x$ (Å)')
ax.set_ylabel('$y$ (Å)')
zlabel = ax.set_zlabel('$z$ (Å)')
zlabel.set_rotation(90)  # Try 90 or 270 for vertical alignment
ax.zaxis.set_rotate_label(False)  # Disable auto-rotation
zlabel.set_verticalalignment('bottom')  # Adjust vertical position
zlabel.set_horizontalalignment('center')

ax.view_init(elev=30, azim=225)

plt.tight_layout()
plt.savefig('monoIII-300.pdf', bbox_inches='tight', pad_inches=0.35)
#plt.show()
plt.close()
##########################

"""
####### MONOIII 1200K #######
fig = plt.figure(figsize=(4, 4))
ax = fig.add_subplot(111, projection='3d')

ax.set_title('$T=1200$ K')

path = 'data/monoIII/1200K/XDATCAR'
number_steps = 30000
nblocks = 10
time_step = 0.0015

atom_subcell = 0 # to denote the subcell in the supercell

atoms_Hf = []
atom = 1 + atom_subcell
for _ in range(8):
    atoms_Hf.append(atom)
    atom = atom + 18

atoms_O = []
atom = 145 + atom_subcell
for _ in range(16):
    atoms_O.append(atom)
    atom = atom + 18

for atom in atoms_Hf:
    plot_traj(path, atom, number_steps, nblocks, time_step, ['brown', 'green'])

for atom in atoms_O:
    plot_traj(path, atom, number_steps, nblocks, time_step, ['red', 'yellow'])

ax.set_xlabel('$x$ (Å)')
ax.set_ylabel('$y$ (Å)')
zlabel = ax.set_zlabel('$z$ (Å)')
zlabel.set_rotation(90)  # Try 90 or 270 for vertical alignment
ax.zaxis.set_rotate_label(False)  # Disable auto-rotation
zlabel.set_verticalalignment('bottom')  # Adjust vertical position
zlabel.set_horizontalalignment('center')

ax.view_init(elev=30, azim=225)

plt.tight_layout()
plt.savefig('monoIII-1200.pdf', bbox_inches='tight', pad_inches=0.35)
#plt.show()
plt.close()
##########################


##### CUBIC diffusing ####
fig = plt.figure(figsize=(4, 4))
ax = fig.add_subplot(111, projection='3d')

path = 'data/monoIII/300K/XDATCAR'
number_steps = 30000
nblocks = 10
time_step = 0.0015

atom_subcell = 0 # to denote the subcell in the supercell

atoms_O = []
atom = 145 + atom_subcell
for _ in range(16):
    atoms_O.append(atom)
    atom = atom + 18

plot_traj(path, atoms_O[0], number_steps, nblocks, time_step, ['red', 'yellow'])

ax.view_init(elev=30, azim=225)

plt.tight_layout()
plt.savefig('monoIII-300-O.pdf', bbox_inches='tight', pad_inches=0.35)
#plt.show()
plt.close()

fig = plt.figure(figsize=(4, 4))
ax = fig.add_subplot(111, projection='3d')

path = 'data/monoIII/1200K/XDATCAR'
number_steps = 30000
nblocks = 10
time_step = 0.0015

atom_subcell = 0 # to denote the subcell in the supercell

atoms_O = []
atom = 145 + atom_subcell
for _ in range(16):
    atoms_O.append(atom)
    atom = atom + 18

plot_traj(path, atoms_O[0], number_steps, nblocks, time_step, ['red', 'yellow'])

ax.view_init(elev=30, azim=225)

plt.tight_layout()
plt.savefig('monoIII-1200-O.pdf', bbox_inches='tight', pad_inches=0.35)
#plt.show()
plt.close()




fig = plt.figure(figsize=(4, 4))
ax = fig.add_subplot(111, projection='3d')

path = 'data/monoIII/300K/XDATCAR'
number_steps = 30000
nblocks = 10
time_step = 0.0015

atom_subcell = 0 # to denote the subcell in the supercell

atoms_Hf = []
atom = 1 + atom_subcell
for _ in range(8):
    atoms_Hf.append(atom)
    atom = atom + 18

plot_traj(path, atoms_Hf[2], number_steps, nblocks, time_step, ['brown', 'green'])

ax.view_init(elev=30, azim=225)

plt.tight_layout()
plt.savefig('monoIII-300-Hf.pdf', bbox_inches='tight', pad_inches=0.35)
#plt.show()
plt.close()

fig = plt.figure(figsize=(4, 4))
ax = fig.add_subplot(111, projection='3d')

path = 'data/monoIII/1200K/XDATCAR'
number_steps = 30000
nblocks = 10
time_step = 0.0015

atom_subcell = 0 # to denote the subcell in the supercell

atoms_Hf = []
atom = 1 + atom_subcell
for _ in range(8):
    atoms_Hf.append(atom)
    atom = atom + 18

plot_traj(path, atoms_Hf[2], number_steps, nblocks, time_step, ['brown', 'green'])

ax.view_init(elev=30, azim=225)

plt.tight_layout()
plt.savefig('monoIII-1200-Hf.pdf', bbox_inches='tight', pad_inches=0.35)
#plt.show()
plt.close()
##########################
"""
