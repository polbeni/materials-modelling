# Pol Benítez Colominas, March 2025
# Universitat Politècnica de Catalunya

# Reads VASP data and creates a file xyz to train MACE

import glob
import random

import numpy as np

from pymatgen.io.vasp.outputs import Vasprun


def create_xyz(vaspruns, prop_train, n_blocks):
    """
    Creates a train and test xyz files to retrain MACE (with train and test randomly splitted)

    Inputs:
    vaspruns: array with all the paths to vasprun.xml files
    prop_train: proportion of training set (proportion for test set will be 1-prop_train)
    n_blocks: save every n_blocks step from AIMD
    """

    # Create variables to store information
    num_atoms = []
    lattice_array = []
    energy_array = []
    stress_array = []
    atoms_array = []
    positions_array = []
    forces_array = []

    # Read for all the vasprun files
    for vasp_file in vaspruns:
        print(f'Processing: {vasp_file}')

        # Try to open a given vasprun.xml file
        try:
            # Try to load those unfinished relaxations as well
            vasprun = Vasprun(vasp_file, exception_on_bad_xml=False)
        except:
            print('Error: vasprun not correctly loaded.')
            continue

        # If there is no problem extract the information
        # Check if it is spc or AIMD (if AIMD apply n_blocks condition)
        for ionic_step_idx, ionic_step in enumerate(vasprun.ionic_steps):
            if (ionic_step_idx%n_blocks == 0) or (ionic_step_idx == 0):
                # Extract data from each ionic step
                structure = ionic_step['structure']
                energy    = ionic_step['e_fr_energy']
                forces    = ionic_step['forces']
                stress    = ionic_step['stress']

                # Change stress units to eV/A^3 (by default VASP returns kBar, see https://github.com/ACEsuit/mace/discussions/542)
                stress = np.array(stress)
                stress = stress / 1602

                # Save all the data in the arrays
                num_atoms.append(len(structure))

                lattice = ' '.join(map(str, structure.lattice.matrix.flatten()))
                lattice_array.append(lattice)

                energy_array.append(energy)

                stress_str = ' '.join(map(str, stress.flatten()))
                stress_array.append(stress_str)

                atoms_step = []
                positions_step = []
                forces_step = []
                for idx, _ in enumerate(structure):
                    atoms_step.append(str(structure[idx].specie))
                    positions_step.append(" ".join(map(str, structure[idx].coords)))
                    forces_step.append(" ".join(map(str, forces[idx])))

                atoms_array.append(atoms_step)
                positions_array.append(positions_step)
                forces_array.append(forces_step)

    # Save all the data in the xyz file
    num_cases = len(num_atoms)

    num_cases_train = int(num_cases * prop_train)

    disorded_struc = create_shuffled_array(num_cases)

    train_file = open('train.xyz', 'w')
    test_file = open('test.xyz', 'w')

    it_structure = 0
    for structure in disorded_struc:
        if it_structure <= (num_cases_train - 1):
            # Write the atoms number
            train_file.write(f'{num_atoms[structure]}\n')

            # Write the metadata
            train_file.write(f'Lattice=\"{lattice_array[structure]}\" Properties=species:S:1:pos:R:3:forces:R:3 energy={energy_array[structure]} REF_stress=\"{stress_array[structure]}\"\n')

            # Write the atom data
            for num_atom in range(num_atoms[structure]):
                train_file.write(f'{atoms_array[structure][num_atom]} {positions_array[structure][num_atom]} {forces_array[structure][num_atom]}\n')
        else:
            # Write the atoms number
            test_file.write(f'{num_atoms[structure]}\n')

            # Write the metadata
            test_file.write(f'Lattice=\"{lattice_array[structure]}\" Properties=species:S:1:pos:R:3:forces:R:3 energy={energy_array[structure]} REF_stress=\"{stress_array[structure]}\"\n')

            # Write the atom data
            for num_atom in range(num_atoms[structure]):
                test_file.write(f'{atoms_array[structure][num_atom]} {positions_array[structure][num_atom]} {forces_array[structure][num_atom]}\n')

        it_structure = it_structure + 1

    train_file.close()
    test_file.close()


def create_shuffled_array(n):
    """
    Returns a list with the intergers from 0 to n-1 disordered

    Inputs:
        n: size of the list
    """

    numbers = list(range(n))
 
    random.shuffle(numbers)

    return numbers

### INPUT PARAMETERS ###
n_blocks = 5 # Save every n_blocks step from AIMD
train_ratio = 0.8 # Proportion of data for training
path_to_vaspruns = 'data/**/vasprun.xml'

# Take the vasprun.xml data from the specified path
vasp_files = glob.glob(path_to_vaspruns, recursive=True)  # Modify path

# Create the xyz files
create_xyz(vasp_files, train_ratio, n_blocks)