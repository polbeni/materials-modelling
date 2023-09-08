def generate_POSCAR(input_file, num_chemical, num_config, output_file):
    """
    This function generates a POSCAR from the XDATCAR file
    
    input_file: name or path of the input file, it should be a XDATCAR file
    num_chemical: number of different chemical species in the POSCAR
    num_config: number of the configuration at the XDATCAR file
    output_file: name of the output file, a new POSCAR with the positions of the num_config configuration
    """

    XDATCAR = open(input_file, "r")
    
    for x in range(7):
        number_atoms_line = XDATCAR.readline()

    number_atoms = 0
    for x in range(num_chemical):
        number_atoms = number_atoms + int(number_atoms_line.split()[x])

    XDATCAR.close()

    XDATCAR = open(input_file, "r")
    new_POSCAR = open(output_file, "w")

    for x in range(7):
        XDATCAR_line = XDATCAR.readline()
        if x < 5:
            new_POSCAR.writelines(XDATCAR_line)
        else:
            new_POSCAR.writelines(str(XDATCAR_line.split()[0]) + ' ' + str(XDATCAR_line.split()[1]) + ' ' + str(XDATCAR_line.split()[2]) + '\n')

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


num_POSCAR = 30000

for x in range(10):
    name_file = 'generated_POSCAR/POSCAR-' + "{:03d}".format(x+1)
    
    num_POSCAR = num_POSCAR + (60000 - 30000)/10

    print('POSCAR-' + "{:03d}".format(x+1) + ':')
    generate_POSCAR('XDATCAR_merged', 3, num_POSCAR, name_file)
    print(' ')
