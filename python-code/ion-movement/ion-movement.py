import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib.ticker import MultipleLocator

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

    for x in range(7):
        actual_line = XDATCAR_file.readline()
    
    num_atoms = int(actual_line.split()[0]) + int(actual_line.split()[1]) + int(actual_line.split()[2])

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
                ion_movement[num_configuration, 1] = float(actual_line.split()[0]) * 19.17220921 # lattice parameter
                ion_movement[num_configuration, 2] = float(actual_line.split()[1]) * 19.17220921
                ion_movement[num_configuration, 3] = float(actual_line.split()[2]) * 19.17220921

                num_configuration = num_configuration + 1
        

        

        

    XDATCAR_file.close()

    return ion_movement

# 1 ion
ion_num = 11
ion_m = get_movement('XDATCAR-600', ion_num, 25000, 10, 0.0015)

ion_m[:,1] = ion_m[:,1] - ion_m[0,1]
ion_m[:,2] = ion_m[:,2] - ion_m[0,2]
ion_m[:,3] = ion_m[:,3] - ion_m[0,3]

min_x = np.min(ion_m[:,1]) - 0.15*(np.max(ion_m[:,1]) - np.min(ion_m[:,1]))
max_x = np.max(ion_m[:,1]) + 0.15*(np.max(ion_m[:,1]) - np.min(ion_m[:,1]))

min_y = np.min(ion_m[:,2]) - 0.15*(np.max(ion_m[:,2]) - np.min(ion_m[:,2]))
max_y = np.max(ion_m[:,2]) + 0.15*(np.max(ion_m[:,2]) - np.min(ion_m[:,2]))

min_z = np.min(ion_m[:,3]) - 0.15*(np.max(ion_m[:,3]) - np.min(ion_m[:,3]))
max_z = np.max(ion_m[:,3]) + 0.15*(np.max(ion_m[:,3]) - np.min(ion_m[:,3]))

fig = plt.figure()
ax = fig.add_subplot(projection='3d')

ax.set_xlabel('$\\Delta$x (Å)', fontsize=20, labelpad=10)
ax.set_ylabel('$\\Delta$y (Å)', fontsize=20, labelpad=10)
ax.set_zlabel('')  # Remove the default z-axis label

# Manually add the z-axis label with rotation
ax.text2D(-0.1, 0.5, '$\\Delta$z (Å)', fontsize=20, rotation=90, va='center', ha='center', transform=ax.transAxes)

ax.text2D(0.7, 0.75, 'Ag$_3$SBr', fontsize=22, va='center', ha='center', transform=ax.transAxes)
ax.text2D(0.7, 0.65, 'T=600 K', fontsize=22, va='center', ha='center', transform=ax.transAxes)


#ax.set_xlim((-0.15, 0.15))
#ax.set_ylim((-0.15, 0.15))
#ax.set_zlim((-0.15, 0.15))

ax.set_xticklabels(ion_m[:,1], fontsize=17)
ax.set_yticklabels(ion_m[:,2], fontsize=17)
ax.set_zticklabels(ion_m[:,3], fontsize=17)

formatter = FuncFormatter(lambda x, _: "{:g}".format(x))
ax.xaxis.set_major_formatter(formatter)
ax.yaxis.set_major_formatter(formatter)
ax.zaxis.set_major_formatter(formatter)


scatter = ax.scatter(ion_m[:,1], ion_m[:,2], ion_m[:,3], c=ion_m[:,0])

major_locator = MultipleLocator(1)  
ax.xaxis.set_major_locator(major_locator)
ax.yaxis.set_major_locator(major_locator)
ax.zaxis.set_major_locator(major_locator)

# Fixing camera position
ax.view_init(elev=15, azim=35)  # Set desired elevation and azimuth

# Rotate z-axis labels
ax.zaxis.label.set_rotation(90)

# Adjust tick label padding
ax.tick_params(axis='x', pad=0)
ax.tick_params(axis='y', pad=0)
ax.tick_params(axis='z', pad=4)


plt.savefig('Ag3SBr-600.pdf')

plt.tight_layout()
#plt.show()




#ion_num = 11
ion_m = get_movement('XDATCAR-200', ion_num, 25000, 10, 0.0015)

ion_m[:,1] = ion_m[:,1] - ion_m[0,1]
ion_m[:,2] = ion_m[:,2] - ion_m[0,2]
ion_m[:,3] = ion_m[:,3] - ion_m[0,3]

fig = plt.figure()
ax = fig.add_subplot(projection='3d')

ax.set_xlabel('$\\Delta$x (Å)', fontsize=20, labelpad=10)
ax.set_ylabel('$\\Delta$y (Å)', fontsize=20, labelpad=10)
ax.set_zlabel('')  # Remove the default z-axis label

# Manually add the z-axis label with rotation
ax.text2D(-0.15, 0.42, '$\\Delta$z (Å)', fontsize=20, rotation=90, va='center', ha='center', transform=ax.transAxes)

ax.text2D(0.7, 0.8, 'Ag$_3$SBr', fontsize=22, va='center', ha='center', transform=ax.transAxes)
ax.text2D(0.7, 0.7, 'T=200 K', fontsize=22, va='center', ha='center', transform=ax.transAxes)

#ax.set_xlim((-0.07, 0.07))
#ax.set_ylim((-0.07, 0.07))
#ax.set_zlim((-0.07, 0.07))

ax.set_xticklabels(ion_m[:,1], fontsize=17)
ax.set_yticklabels(ion_m[:,2], fontsize=17)
ax.set_zticklabels(ion_m[:,3], fontsize=17)

formatter = FuncFormatter(lambda x, _: "{:g}".format(x))
ax.xaxis.set_major_formatter(formatter)
ax.yaxis.set_major_formatter(formatter)
ax.zaxis.set_major_formatter(formatter)

scatter = ax.scatter(ion_m[:,1], ion_m[:,2], ion_m[:,3], c=ion_m[:,0])

major_locator = MultipleLocator(1)  
ax.xaxis.set_major_locator(major_locator)
ax.yaxis.set_major_locator(major_locator) 
ax.zaxis.set_major_locator(major_locator)

# Fixing camera position
ax.view_init(elev=45, azim=35)  # Set desired elevation and azimuth

# Rotate z-axis labels
ax.zaxis.label.set_rotation(90)

# Adjust tick label padding
ax.tick_params(axis='x', pad=0)
ax.tick_params(axis='y', pad=0)
ax.tick_params(axis='z', pad=8)

plt.savefig('Ag3SBr-200.pdf')

plt.tight_layout()
plt.show()



fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.set_xlabel('$\\Delta$x (Å)')
ax.set_ylabel('$\\Delta$y (Å)')
ax.set_zlabel('$\\Delta$z (Å)')
#ax.set_xlim((min_x, max_x))
#ax.set_ylim((min_y, max_y))
#ax.set_zlim((min_z, max_z))
formatter = FuncFormatter(lambda x, _: "{:g}".format(x))
ax.xaxis.set_major_formatter(formatter)
ax.yaxis.set_major_formatter(formatter)
ax.zaxis.set_major_formatter(formatter)
#ax.set_xlim(0, 1)
#ax.set_ylim(0, 1)
#ax.set_zlim(0, 1)
scatter = ax.scatter(ion_m[:,1], ion_m[:,2], ion_m[:,3], c=ion_m[:,0])
cbar = plt.colorbar(scatter)
cbar.ax.tick_params(labelsize=15) 
cbar.set_label('t (ps)', fontsize=15)
#plt.savefig('ion_movement.pdf')

plt.tight_layout()
#plt.show()
