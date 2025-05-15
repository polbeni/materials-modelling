# Pol Benítez Colominas, May 2025
# Universitat Politècnica de Catalunya

# Study the movement and rotation of molecules in a molecular crystal

import numpy as np
import matplotlib.pyplot as plt

def get_num_atoms(path_XDATCAR):
    """
    Gets the total number of atoms and carbons in the AIMD simulation

    Inputs:
        path_XDATCAR: path to the XDATCAR file
    """

    XDATCAR = open(path_XDATCAR, 'r')

    for _ in range(6):
        XDATCAR.readline()

    line_num_atoms = XDATCAR.readline()

    num_atoms = 0
    for x in range(len(line_num_atoms.split())):
        num_atoms = num_atoms + int(line_num_atoms.split()[x])
    
    num_carbons = int(line_num_atoms.split()[0])

    XDATCAR.close()

    return num_atoms, num_carbons

def wrap_coordinates(x, y, z):
    """
    Ensure periodic boundary conditions

    Inputs:
        x, y, z: coordinates
    """

    return [x % 1, y % 1, z % 1]

def generate_pairs(total_num_carbons):
    """
    Generates the different pair of atom indices representing each molecule

    Inputs:
        total_num_carbons: total number of carbons in the AIMD simulation
    """

    pairs_arr = []

    for number in range(int(total_num_carbons/4)):
        pair = [number + 1, int(total_num_carbons/2) + number + 1]
        pairs_arr.append(pair)

    return pairs_arr

def generate_molecules(total_num_atoms):
    """
    Generates the different atom indices representing each molecule

    Inputs:
        total_num_atoms: total number of atoms in the AIMD simulation
    """

    positions_list = []

    number_atoms_molecule = 10

    for number in range(int(total_num_atoms/number_atoms_molecule)):
        positions = [number + 1, number + int(total_num_atoms/number_atoms_molecule)*1 + 1, number + int(total_num_atoms/number_atoms_molecule)*2 + 1, number + int(total_num_atoms/number_atoms_molecule)*3 + 1,
                number + int(total_num_atoms/number_atoms_molecule)*4 + 1, number + int(total_num_atoms/number_atoms_molecule)*5 + 1, number + int(total_num_atoms/number_atoms_molecule)*6 + 1, number + int(total_num_atoms/number_atoms_molecule)*7 + 1,
                number + int(total_num_atoms/number_atoms_molecule)*8 + 1, number + int(total_num_atoms/number_atoms_molecule)*9 + 1]
        positions_list.append(positions)

    return positions_list

def center_mass_position(positions, masses):
    """
    Determines the position of the center of mass of a molecule

    Inputs:
        positions: array with the atoms positions
        masses: array with the atoms masses
    """

    total_mass = 0
    for mass in masses:
        total_mass = total_mass + mass

    x_cm = 0
    y_cm = 0
    z_cm = 0
    for atom in range(len(positions)):
        x_cm = x_cm + (masses[atom] * positions[atom][0])
        y_cm = y_cm + (masses[atom] * positions[atom][1])
        z_cm = z_cm + (masses[atom] * positions[atom][2])

    vec_cm = [x_cm / total_mass, y_cm / total_mass, z_cm / total_mass]

    return vec_cm

def get_vectors_AIMD(path_XDATCAR, num_atoms, num_iterations, atom_1, atom_2, molecule_positions, masses_list, lattice_parameters):
    """
    Returns a list with the vector for a pair of atoms and the position of the center of mass for each time step of the AIMD

    Inputs:
        path_XDATCAR: path to the XDATCAR file
        num_atoms: total number of atoms in the AIMD simulation
        num_iterations: total number of iterations to include
        atom_1, atom_2: indices of the two considered atoms
        masses_list: array with all the masses in the molecule
        lattice_parameters: lattice parameters of the unit cell
    """

    XDATCAR = open(path_XDATCAR, 'r')

    for _ in range(7):
        XDATCAR.readline()

    list_vectors = []
    list_positions = []

    reference_position = [0, 0, 0]

    for it_md in range(num_iterations):
        XDATCAR.readline()

        positions_atoms = []
        for it_atom in range(num_atoms):
            line = XDATCAR.readline()
            
            if (it_atom + 1) == atom_1:
                vec1 = [float(line.split()[0])*lattice_parameters[0], float(line.split()[1])*lattice_parameters[1], float(line.split()[2])*lattice_parameters[2]]
            elif (it_atom + 1) == atom_2:
                vec2 = [float(line.split()[0])*lattice_parameters[0], float(line.split()[1])*lattice_parameters[1], float(line.split()[2])*lattice_parameters[2]]

            if (it_atom + 1) in molecule_positions:
                positions_atoms.append([float(line.split()[0])*lattice_parameters[0], float(line.split()[1])*lattice_parameters[1], float(line.split()[2])*lattice_parameters[2]])
        
        if it_md == 0:
            reference_position = define_origin_molecule(positions_atoms)
        else:
            reference_position = cm_point

        corrected_positions = correct_positions(positions_atoms, reference_position, lattice_parameters)
        vec1 = corrected_positions[0]
        vec2 = corrected_positions[1]
        final_vector = [vec2[0] - vec1[0], vec2[1] - vec1[1], vec2[2] - vec1[2]]
        list_vectors.append(final_vector)

        corrected_positions = correct_positions(positions_atoms, reference_position, lattice_parameters)
        cm_point = center_mass_position(corrected_positions, masses_list)
        list_positions.append(cm_point)

    return list_vectors, list_positions

def position_change(list_positions):
    """
    Computes the change of position with respect the initial position

    Inputs:
        list_positions: array with the positions
    """

    inital_position_x = list_positions[0][0]
    inital_position_y = list_positions[0][1]
    inital_position_z = list_positions[0][2]

    disp_x = []
    disp_y = []
    disp_z = []

    for it_num in range(len(list_positions)):
        disp_x.append(list_positions[it_num][0] - inital_position_x)
        disp_y.append(list_positions[it_num][1] - inital_position_y)
        disp_z.append(list_positions[it_num][2] - inital_position_z)

    return disp_x, disp_y, disp_z

def compute_angle(vec1, vec2):
    """
    Determines the angle between two vectors (in deg)

    Inputs:
        vec1, vec2: the two vectors
    """

    angle = np.arccos((vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2]) / (np.sqrt((vec1[0]**2) + (vec1[1]**2) + (vec1[2]**2))*np.sqrt((vec2[0]**2) + (vec2[1]**2) + (vec2[2]**2))))

    return np.degrees(angle)

def compute_angle_360(vec1, vec2):
    """
    Determines the angle between two vectors (in deg) with the total 360 deg

    Inputs:
        vec1, vec2: the two vectors (2d)
    """

    vec1 = np.array(vec1)
    vec2 = np.array(vec2)
    
    # Normalize
    v1_u = vec1 / np.linalg.norm(vec1)
    v2_u = vec2 / np.linalg.norm(vec2)
    
    # Dot and cross
    dot = np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)
    cross = v1_u[0]*v2_u[1] - v1_u[1]*v2_u[0]  # equivalent to z-component of 2D cross

    angle_rad = np.arccos(dot)
    if cross < 0:
        angle_rad = 2 * np.pi - angle_rad

    return np.degrees(angle_rad)

def angular_change(list_vectors):
    """
    Computes the change of the angular coordinates

    Inputs:
        list_vectors: array with the vectors
    """

    theta = []
    phi = []

    for it_num in range(len(list_vectors)):
        angle = compute_angle_360([list_vectors[it_num][0], list_vectors[it_num][1]], [1, 0])
        theta.append(angle)

        angle = compute_angle([list_vectors[it_num][0], list_vectors[it_num][1], list_vectors[it_num][2]], [0, 0, 1])
        phi.append(angle)

    return theta, phi

def make_plots(time_step, num_steps, molecule_disp_x, molecule_disp_y, molecule_disp_z, molecule_rot_theta, molecule_rot_phi, type_struc, temp):
    """
    Generate a displacements and rotations plot with the results

    Inputs:
        time_step: time step used in the AIMD
        num_steps: total number of steps in the AIMD
        molecule_...: arrays with displacements and rotations
        type_struc: type of the structure: vacancy, pristine...
        temp: temperature of the simulation
    """
    steps = np.linspace(0, time_step*num_steps, num_steps)

    # Plot displacements
    plt.figure()
    fig, axs = plt.subplots(3, 1, figsize=(4, 5))

    axs[0].set_title(f'{type_struc} T={temp} K')

    axs[0].set_ylabel('$\\Delta x$ ($\\AA$)')
    axs[1].set_ylabel('$\\Delta y$ ($\\AA$)')
    axs[2].set_ylabel('$\\Delta z$ ($\\AA$)')
    axs[2].set_xlabel('time (ps)')
    
    for molecule in range(len(molecule_disp_x)):
        axs[0].plot(steps, molecule_disp_x[molecule])
        axs[1].plot(steps, molecule_disp_y[molecule])
        axs[2].plot(steps, molecule_disp_z[molecule])

    plt.tight_layout()
    plt.savefig(f'resultats/displacements-{type_struc}-{temp}.pdf')

    # Plot rotations
    plt.figure()
    fig, axs = plt.subplots(2, 1, figsize=(4, 4))

    axs[0].set_title(f'{type_struc} T={temp} K')

    axs[0].set_ylabel('$\\theta$')
    axs[1].set_ylabel('$\\phi$')
    axs[1].set_xlabel('time (ps)')

    steps = np.linspace(0, time_step*num_steps, num_steps)
    for molecule in range(len(molecule_rot_theta)):
        angles_deg = molecule_rot_theta[molecule]
        angles_rad = np.radians(angles_deg)
        unwrapped_rad = np.unwrap(angles_rad)
        unwrapped_deg = np.degrees(unwrapped_rad)
        axs[0].plot(steps, unwrapped_deg)

        angles_deg = molecule_rot_phi[molecule]
        angles_rad = np.radians(angles_deg)
        unwrapped_rad = np.unwrap(angles_rad)
        unwrapped_deg = np.degrees(unwrapped_rad)
        axs[1].plot(steps, unwrapped_deg)

    axs[0].set_ylim(0, 360)
    axs[1].set_ylim(0, 180)

    plt.tight_layout()
    plt.savefig(f'resultats/rotations-{type_struc}-{temp}.pdf')

def define_origin_molecule(atoms_positions):
    """
    It defines an original positions for all the molecules

    Inputs:
        atoms_positions: list with the positions of all the atoms
    """

    center = atoms_positions[0]

    return center

def correct_positions(atoms_positions, center_molecule, lattice_parameters, threshold=8):
    """
    Corrects the positions of the atoms that are in the other side of the pbc of the cell

    Inputs:
        atoms_positions: list with the positions of all the atoms
        center_molecule: the center of the molecule in the previous step, it is used as a reference
        lattice_parameters: lattice parameters of the unit cell
        threshold: distance threshold to see if the atom is not in the correct side
    """
    new_positions = []

    for atom in atoms_positions:
        new_pos = []

        if (atom[0] - center_molecule[0]) > threshold:
            new_pos.append(atom[0] - lattice_parameters[0])
        elif abs(-atom[0] + center_molecule[0]) > threshold:
            new_pos.append(atom[0] + lattice_parameters[0])
        else:
            new_pos.append(atom[0])

        if (atom[1] - center_molecule[1]) > threshold:
            new_pos.append(atom[1] - lattice_parameters[1])
        elif abs(-atom[1] + center_molecule[1]) > threshold:
            new_pos.append(atom[1] + lattice_parameters[1])
        else:
            new_pos.append(atom[1])

        if (atom[2] - center_molecule[2]) > threshold:
            new_pos.append(atom[2] - lattice_parameters[2])
        elif abs(-atom[2] + center_molecule[2]) > threshold:
            new_pos.append(atom[2] + lattice_parameters[2])
        else:
            new_pos.append(atom[2])

        new_positions.append(new_pos)

    return new_positions

def run_everything(path, structure, temperature, lattice_parameters):
    """
    It runs all the steps

    Inputs:
        path: path to the XDATCAR file
        structure: type of the structure: vacancy, pristine...
        temperature: temperature of the simulation
        lattice_parameters: lattice parameters of the unit cell
    """

    # Determine total number of atoms and generate the carbon pairs
    num_atoms, num_carbons = get_num_atoms(path)
    pairs_arr = generate_pairs(num_carbons)
    molecules_arr = generate_molecules(num_atoms)

    # Save the vectors and positions for all the pairs
    pair_vectors = []
    pair_positions = []

    masses_list = [12.011, 12.011, 12.011, 12.011, 1.008, 1.008, 1.008, 1.008, 14.007, 14.007]

    it_molecule = 0
    for pair in pairs_arr:
        print(f'Saving pair {pair}')
        vectors, positions = get_vectors_AIMD(path, num_atoms, num_steps, pair[0], pair[1], molecules_arr[it_molecule], masses_list, lattice_parameters)

        pair_vectors.append(vectors)
        pair_positions.append(positions)

        it_molecule = it_molecule + 1

    # Compute the displacements and angular changes for each molecule
    molecule_disp_x = []
    molecule_disp_y = []
    molecule_disp_z = []
    molecule_rot_theta = []
    molecule_rot_phi = []

    for molecule in range(len(pair_vectors)):
        print(f'Computing molecule {molecule + 1}')
        disp_x, disp_y, disp_z = position_change(pair_positions[molecule])
        theta, phi = angular_change(pair_vectors[molecule])

        molecule_disp_x.append(disp_x)
        molecule_disp_y.append(disp_y)
        molecule_disp_z.append(disp_z)
        molecule_rot_theta.append(theta)
        molecule_rot_phi.append(phi)

    make_plots(time_step, num_steps, molecule_disp_x, molecule_disp_y, molecule_disp_z, molecule_rot_theta, molecule_rot_phi, structure, temperature)


# Global parameters
num_steps = 20000
time_step = 0.005
lattice_parameters = [18.182, 16.553, 16.931]

# Run simulation
path = 'total-XDATCAR-files/pristine/T-200/XDATCAR'
structure = 'Pristine'
temperature = '200'

run_everything(path, structure, temperature, lattice_parameters)

path = 'total-XDATCAR-files/pristine/T-400/XDATCAR'
structure = 'Pristine'
temperature = '400'

run_everything(path, structure, temperature, lattice_parameters)

path = 'total-XDATCAR-files/vacancy-1/T-200/XDATCAR'
structure = 'Vacancy-1'
temperature = '200'

run_everything(path, structure, temperature, lattice_parameters)

path = 'total-XDATCAR-files/vacancy-1/T-400/XDATCAR'
structure = 'Vacancy-1'
temperature = '400'

run_everything(path, structure, temperature, lattice_parameters)

path = 'total-XDATCAR-files/vacancy-2/T-200/XDATCAR'
structure = 'Vacancy-2'
temperature = '200'

run_everything(path, structure, temperature, lattice_parameters)

path = 'total-XDATCAR-files/vacancy-2/T-400/XDATCAR'
structure = 'Vacancy-2'
temperature = '400'

run_everything(path, structure, temperature, lattice_parameters)

path = 'total-XDATCAR-files/vacancy-3/T-200/XDATCAR'
structure = 'Vacancy-3'
temperature = '200'

run_everything(path, structure, temperature, lattice_parameters)

path = 'total-XDATCAR-files/vacancy-3/T-400/XDATCAR'
structure = 'Vacancy-3'
temperature = '400'

run_everything(path, structure, temperature, lattice_parameters)

path = 'total-XDATCAR-files/vacancy-4/T-200/XDATCAR'
structure = 'Vacancy-4'
temperature = '200'

run_everything(path, structure, temperature, lattice_parameters)

path = 'total-XDATCAR-files/vacancy-4/T-400/XDATCAR'
structure = 'Vacancy-4'
temperature = '400'

run_everything(path, structure, temperature, lattice_parameters)