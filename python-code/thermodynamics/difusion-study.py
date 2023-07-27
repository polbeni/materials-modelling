import numpy as np
import matplotlib.pyplot as plt

def generate_positions_for_ion(XDATCARfile, number_chemical, number_steps, number_jump, number_ion, output_file):
    """
    This function generates a file with the positions of an ion during some iterations

    XDATCARfile: the path to the XDATCAR file
    number chemical: number of different chemical species in the POSCAR
    number_steps: the number of steps that we want to take from the XDATCAR file
    number_jump: the number of steps that we want to save, example: if number_steps=1000 and
                 number_jump=10, the function will only storage 100 steps
    number_ion: the number ion in the POSCAR that we want to know  
    output_file: name of the output file
    """
    
    XDATCAR = open(XDATCARfile, "r")
    
    for x in range(7):
        number_atoms_line = XDATCAR.readline()

    number_atoms = 0
    for x in range(number_chemical):
        number_atoms = number_atoms + int(number_atoms_line.split()[x])

    XDATCAR.close()

    XDATCAR = open(XDATCARfile, "r")
    POSITIONS = open(output_file, "w")
    
    finish_condition = False
    num_it = 0

    actual_iteration = 1
    iteration_step = 7 + 1 + number_ion
    while finish_condition == False:
        actual_line = XDATCAR.readline()

        if actual_iteration == iteration_step:
            POSITIONS.writelines(actual_line)

            num_it = num_it + 1
            iteration_step = iteration_step + (number_atoms + 1)*number_jump

        if num_it >= (number_steps/number_jump):
            finish_condition = True

        actual_iteration = actual_iteration + 1

    XDATCAR.close()
    POSITIONS.close()

    return print(f'File with positions generated for the ion number {number_ion}')

total_steps = 60000

generate_positions_for_ion('XDATCAR', 3, total_steps, 100, 4, 'ionPOSITIONS1')
generate_positions_for_ion('XDATCAR', 3, total_steps, 100, 26, 'ionPOSITIONS2')
generate_positions_for_ion('XDATCAR', 3, total_steps, 100, 39, 'ionPOSITIONS3')

data_ion1 = np.genfromtxt('ionPOSITIONS1', delimiter='')
data_ion2 = np.genfromtxt('ionPOSITIONS2', delimiter='')
data_ion3 = np.genfromtxt('ionPOSITIONS3', delimiter='')

time_array = np.linspace(0, total_steps*0.0015, len(data_ion1))

plt.figure()
fig, (ax1, ax2, ax3) = plt.subplots(3)
fig.suptitle('Cubic Pm-3m 2x2x2, T=600K')
ax1.set_xlabel('t (ps)')
ax1.set_ylabel('x')
ax1.plot(time_array, data_ion1[:,0], label='$Ag^{+}$')
ax1.plot(time_array, data_ion2[:,0], label='$S^{-2}$')
ax1.plot(time_array, data_ion3[:,0], label='$Br^{-}$')
ax1.legend()
ax2.set_xlabel('t (ps)')
ax2.set_ylabel('y')
ax2.plot(time_array, data_ion1[:,1], label='$Ag^{+}$')
ax2.plot(time_array, data_ion2[:,1], label='$S^{-2}$')
ax2.plot(time_array, data_ion3[:,1], label='$Br^{-}$')
#ax2.legend()
ax3.set_xlabel('t (ps)')
ax3.set_ylabel('z')
ax3.plot(time_array, data_ion1[:,2], label='$Ag^{+}$')
ax3.plot(time_array, data_ion2[:,2], label='$S^{-2}$')
ax3.plot(time_array, data_ion3[:,2], label='$Br^{-}$')
plt.tight_layout()
plt.savefig('plots/ion-displacment600.pdf')