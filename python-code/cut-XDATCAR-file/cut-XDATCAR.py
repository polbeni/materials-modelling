# Pol Benítez Colominas, December 2024
# Universitat Politècnica de Catalunya

# Generate a XDATCAR from a desired starting and ending point

def generate_XDATCAR(path_input, path_output, start_step, end_step, num_chemical):
    """
    Generate a XDATCAR with the desired startin and ending point from an original XDATCAR file

    Inputs:
        path_input -> path of the intial step of the XDATCAR file
        path_output -> path of the final step of the XDATCAR file
        start_step -> first step (the steps count starts at 0)
        end_step -> first step (the steps count starts at 0)
        num_chemical -> number of different chemical species in the POSCAR
    """

    XDATCAR_ini = open(path_input, 'r')
    XDATCAR_fin = open(path_output, 'w')

    for x in range(7):
        line = XDATCAR_ini.readline()
        XDATCAR_fin.write(line)

    number_atoms = 0
    for x in range(num_chemical):
        number_atoms = number_atoms + int(line.split()[x])

    it_pos = 0
    it_pos_new = 0
    condition_finish = False
    while condition_finish == False:
        if (it_pos >= start_step) and (it_pos <= end_step):
            for num_line in range(number_atoms + 1):
                line = XDATCAR_ini.readline()

                if num_line == 0:
                    XDATCAR_fin.write(f'{line.split()[0]} {line.split()[1]}   {str(it_pos_new + 1)}\n')
                else:
                    XDATCAR_fin.write(line)

            it_pos_new = it_pos_new + 1
        elif it_pos > end_step:
            condition_finish = True
        else:
            for num_line in range(number_atoms + 1):
                line = XDATCAR_ini.readline()

        it_pos = it_pos + 1

    XDATCAR_ini.close()
    XDATCAR_fin.close()

generate_XDATCAR('XDATCAR1', 'XDATCAR-final', 0, 8300, 2)