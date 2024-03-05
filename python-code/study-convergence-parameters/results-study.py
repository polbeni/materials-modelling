# Pol Benítez Colominas, March 2024
# Universitat Politècnica de Catalunya

# return the energy results of the parameter study

def get_energy(oszicar_file):
    """
    This functions gets the ground state for a VASP simulation results

    Inputs:
        atoms_unit_cell -> number of atoms in the unit cell
        oszicar_file -> path direction for the OSZICAR file
    """


    file = open(oszicar_file, "r")

    with file as f:
        for line in f:
            pass
        last_line = line

    ground_state = float(last_line.split()[4])

    file.close()

    return ground_state

possible_kpoints = [4, 5, 6, 7, 8, 9, 10] # desired kpoints grid
possible_encuts = [450, 500, 550, 600, 650, 700, 750, 800] # desired encut values

for x in possible_kpoints:
    for y in possible_encuts:
        name_path = f'{x:01d}x{x:01d}x{x:01d}-{y:01d}'

        energy = get_energy(name_path + '/OSZICAR')

        print(f'The obtained energy for kpoints {x:01d}x{x:01d}x{x:01d} and ENCUT={y:01d} is {energy:.5f} eV')