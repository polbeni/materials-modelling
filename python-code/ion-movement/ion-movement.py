import numpy as np
import matplotlib.pyplot as plt

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
                ion_movement[num_configuration, 1] = float(actual_line.split()[0])
                ion_movement[num_configuration, 2] = float(actual_line.split()[1])
                ion_movement[num_configuration, 3] = float(actual_line.split()[2])

                num_configuration = num_configuration + 1
        

        

        

    XDATCAR_file.close()

    return ion_movement

# 1 ion
ion_num = 16
ion_m = get_movement('XDATCAR', ion_num, 25000, 10, 0.0015)

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.set_xlabel('x (Å)')
ax.set_ylabel('y (Å)')
ax.set_zlabel('z (Å)')
#ax.set_xlim(0, 1)
#ax.set_ylim(0, 1)
#ax.set_zlim(0, 1)
scatter = ax.scatter(ion_m[:,1], ion_m[:,2], ion_m[:,3], c=ion_m[:,0])
cbar = plt.colorbar(scatter)
cbar.set_label('t (ps)')
#plt.savefig('ion_movement.pdf')
plt.show()


# multiple ions
"""
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.set_xlabel('x (Å)')
ax.set_ylabel('y (Å)')
ax.set_zlabel('z (Å)')
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_zlim(0, 1)
for x in range(0, 190, 19):
    ion_m = get_movement('XDATCAR', x+1, 25000, 100, 0.0015)
    scatter = ax.scatter(ion_m[:,1], ion_m[:,2], ion_m[:,3], c=ion_m[:,0])
cbar = plt.colorbar(scatter)
cbar.set_label('t (ps)')
plt.savefig('ion_movement.pdf')
plt.show()
"""