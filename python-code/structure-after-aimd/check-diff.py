# Pol Benítez Colominas, October 2023
# Universitat Politècnica de Catalunya

import numpy as np

num_silver_ions = 192

prop_matrix = np.zeros((num_silver_ions, 3))

results = open('results', "r")
results.readline()

for x in range(num_silver_ions):
    actual_line = results.readline()

    prop_matrix[x,0] = float(actual_line.split()[3])/float(actual_line.split()[2])
    prop_matrix[x,1] = float(actual_line.split()[6])/float(actual_line.split()[5])
    prop_matrix[x,2] = float(actual_line.split()[9])/float(actual_line.split()[8])

results.close()

prop_disp = open('prop_disp', "w")
prop_disp.write("ion_#     x_diff/x_error     z_diff/z_error     z_diff/z_error\n")

for x in range(num_silver_ions):
    new_line = f'{x+1}     {prop_matrix[x,0]:.4f}     {prop_matrix[x,1]:.4f}     {prop_matrix[x,2]:.4f}\n'

    prop_disp.write(new_line)

prop_disp.close()