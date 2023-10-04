# Pol Benítez Colominas, October 2023
# Universitat Politècnica de Catalunya

from re import A
import numpy as np

def get_positions(path_XDATCAR, ion_number, total_steps, partial_steps):
    """
    This function returns a matrix with the positions of a given atom (ion_number) for each
    configuration reached with aimd simulation from the decided configuration (partial_steps)

    Inputs:
        path_XDATCAR: the path of our XDATCAR file
        ion_number: the position in the POSCAR file of the atom we want to study
        total_steps: number of steps that we want to read in the XDATCAR file
        partial_steps: number of step when we want to start storing the positions
    """
    
    matrix_positions = np.zeros((total_steps - partial_steps, 3))

    XDATCAR_file = open(path_XDATCAR, "r")

    for x in range(7):
        actual_line = XDATCAR_file.readline()
    
    num_atoms = int(actual_line.split()[0]) + int(actual_line.split()[1]) + int(actual_line.split()[2])

    for iteration1 in range(partial_steps - 1):
        XDATCAR_file.readline()

        for iteration2 in range(num_atoms):
            XDATCAR_file.readline()
    
    for iteration1 in range(total_steps - partial_steps):
        XDATCAR_file.readline()

        actual_atom = 1
        for iteration2 in range(num_atoms):
            if actual_atom == ion_number:
                actual_line = XDATCAR_file.readline()

                matrix_positions[iteration1, 0] = float(actual_line.split()[0])
                matrix_positions[iteration1, 1] = float(actual_line.split()[1])
                matrix_positions[iteration1, 2] = float(actual_line.split()[2])
            else:
                XDATCAR_file.readline()

            actual_atom = actual_atom + 1

    XDATCAR_file.close()

    return matrix_positions

XDATCAR_file = open('XDATCAR', "r")

for x in range(7):
    actual_line = XDATCAR_file.readline()
    
num_atoms_POSCAR = int(actual_line.split()[0]) + int(actual_line.split()[1]) + int(actual_line.split()[2])

XDATCAR_file.close()

POSCAR_positions = np.zeros((num_atoms_POSCAR, 3))

POSCAR = open('POSCAR', "r")
for x in range(8):
    POSCAR.readline()

for x in range(num_atoms_POSCAR):
    actual_line = POSCAR.readline()

    POSCAR_positions[x,0] = float(actual_line.split()[0])
    POSCAR_positions[x,1] = float(actual_line.split()[1])
    POSCAR_positions[x,2] = float(actual_line.split()[2])

POSCAR.close()



param_partial = 20000
param_total = 26170

inputs_file = open('results', "w")
# ion_#: ion number
# x,y,z_mean: mean value of x,y,z position for the given ion
# x,y,z_error: statistical error of x,y,z position for the given ion
# x,y,z_diff: difference of the mean x,y,z position for the given ion with respect that in the POSCAR
inputs_file.write('ion_#     x_mean     x_error     x_diff     y_mean     y_error     y_diff     z_mean     z_error     z_diff\n')


NEW_POSCAR = open('NEW_POSCAR', "w")
NEW_POSCAR.write('New structure from aimd simulations\n')

POSCAR = open('POSCAR', "r")
POSCAR.readline()
for x in range(7):
    actual_line = POSCAR.readline()
    new_poscar_string = actual_line
    NEW_POSCAR.write(new_poscar_string)

POSCAR.close()

for x in range(num_atoms_POSCAR):
    positions = get_positions('XDATCAR', x+1, param_total, param_partial)

    mean_x = np.mean(positions[:,0])
    error_x = np.std(positions[:,0])
    diff_x = mean_x - POSCAR_positions[x,0]

    mean_y = np.mean(positions[:,1])
    error_y = np.std(positions[:,1])
    diff_y = mean_y - POSCAR_positions[x,1]

    mean_z = np.mean(positions[:,2])
    error_z = np.std(positions[:,2])
    diff_z = mean_z - POSCAR_positions[x,2]

    inputs_string = f'{x+1}     {mean_x:.10f}     {error_x:.10f}     {diff_x:.10f}     {mean_y:.10f}     {error_y:.10f}     {diff_y:.10f}     {mean_z:.10f}     {error_z:.10f}     {diff_z:.10f}\n'
    inputs_file.write(inputs_string)

    new_poscar_string = f' {mean_x:.10f} {mean_y:.10f} {mean_z:.10f}\n'
    NEW_POSCAR.write(new_poscar_string)

    print(f'Ion computed number: {x+1}')

inputs_file.close()
NEW_POSCAR.close()