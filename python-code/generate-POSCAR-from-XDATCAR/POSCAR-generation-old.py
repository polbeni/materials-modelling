# Pol Benítez Colominas, Sep 2023 - Jan 2025
# Universitat Politècnica de Catalunya

# Get POSCAR structures from a XDATCAR file
# The XDATCAR file should be in the same path

############### INPUTS ###############
initial_step = 10                  # initial step to consider in the XDATCAR
final_step = 40                    # last step to consider in the XDATCAR
number_of_POSCAR = 15              # total number of POSCAR to generate, it has to be < (final_step-initial_step)
path_to_save = 'generated_POSCAR'  # dir to save the generated POSCAR files
######################################

def generate_POSCAR(input_file, num_config, output_file):
    """
    This function generates a POSCAR from the XDATCAR file
    
    input_file: name or path of the input file, it should be a XDATCAR file
    num_config: number of the configuration at the XDATCAR file
    output_file: name of the output file, a new POSCAR with the positions of the num_config configuration
    """

    XDATCAR = open(input_file, "r")
    
    for x in range(7):
        number_atoms_line = XDATCAR.readline()

    number_atoms = 0
    for x in range(len(number_atoms_line.split())):
        number_atoms = number_atoms + int(number_atoms_line.split()[x])

    XDATCAR.close()

    XDATCAR = open(input_file, "r")
    new_POSCAR = open(output_file, "w")

    for x in range(7):
        XDATCAR_line = XDATCAR.readline()
        new_POSCAR.writelines(XDATCAR_line)

    new_POSCAR.writelines('Direct \n')

    finish = False
    desired_config = False
    num_lines_new = 0
    num_line_actual = 7

    while finish == False:
        num_line_actual = num_line_actual + 1
        XDATCAR_line = XDATCAR.readline()

        if desired_config == True:
            new_POSCAR.writelines(XDATCAR_line)
            num_lines_new = num_lines_new + 1

        if (XDATCAR_line.split()[0] == 'Direct'): 
            if (int(XDATCAR_line.split()[2]) == num_config):
                desired_config = True

        if num_lines_new >= number_atoms:
            finish = True

    XDATCAR.close()
    new_POSCAR.close()

    return print(f'POSCAR generated from Direct configuration={num_config} in the provided XDATCAR file')


num_POSCAR = initial_step

for configuration in range(number_of_POSCAR):
    name_file = path_to_save + '/POSCAR-' + "{:05d}".format(configuration + 1)
    
    num_POSCAR = num_POSCAR + int((final_step - initial_step) / number_of_POSCAR)

    print('POSCAR-' + "{:05d}".format(configuration + 1) + ':')
    generate_POSCAR('XDATCAR', num_POSCAR, name_file)
    print(' ')