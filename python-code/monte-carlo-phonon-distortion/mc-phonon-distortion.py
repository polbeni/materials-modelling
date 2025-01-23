# Pol Benítez Colominas, January 2025
# Universitat Politècnica de Catalunya

# Apply harmonic phonon dispersion at given T to a structure

import yaml
import numpy as np


### Physical constants
reduced_planck_ct = 6.5821e-16 # eV·s
boltzmann_ct = 8.6173e-5 # eV·K^-1


### Functions
def bose_einstein(temperature, phonon_energy):
    """
    Computes the Bose-Einstein distribution of a collection of phonons with a given energy

    Inputs:
        temperature -> temperature (in K)
        phonon_energy -> energy of the given phonon (in s^-1 (Hz))
    """

    distribution = 1 / (np.exp((reduced_planck_ct * 2 * np.pi * phonon_energy) / (boltzmann_ct * temperature)) - 1) 

    return distribution

def phonon_amplitude(temperature, phonon_energy, atom_mass):
    """
    Computes the phonon amplitude for a given temperature and phonon mode mu at q

    Inputs:
        temperature -> temperature (in K)
        phonon_energy -> energy of the given phonon (in s^-1 (Hz))
        atom_mass -> atomic mass of the given atom
    """

    dist = bose_einstein(temperature, phonon_energy)

    h_si = reduced_planck_ct * 6.24e-18 # in J·s
    mass = atom_mass * 1.660539e-27 # in Kg

    amplitude = np.sqrt(h_si / (2 * mass * 2 * np.pi * phonon_energy)) * np.sqrt(dist + (1 / 2))

    amplitude = amplitude * 1e10 # in angstrom

    return amplitude

def random_uniform(a, b):
    """
    Generates random uniform number between a and b
    """

    return np.random.uniform(a, b)

def disp_vect(number_atoms, number_q, number_mu, eigenvalues_list, eigenvectors_list, mass_list, temperature):
    """
    Generates a displacement vector for a given ensemble of phonons and at a given temperature

    Inputs:
        number_atoms -> number of atoms in the unit cell
        number_q -> number of points (mesh) in the irreducible Brillouin zone
        number_mu -> number of different phonon modes (3*N, where N is the number of atoms in the unit cell)
        eigenvalues_list -> phonon eigenvalues with format [[mu_1, mu_2, ...]_q_1, [mu_1, mu_2, ...]_q_2, ...]
        eigenvectors_list -> phonon eigenvectors with format [[[atom_1, atom_2, ..]mu_1, [...]mu_2, ...]_q_1, [mu_1[...], mu_2[...], ...]_q_2, ...]
        mass_list -> mass of all the atoms considered
        temperature -> temperature (in K)
    """

    vector = []
    for _ in range(number_atoms):
        vector.append([0, 0, 0])

    for q_val in range(number_q):
        for mu_val in range(number_mu):
            eigenvalue = eigenvalues_list[q_val][mu_val]
            eigenvectors =  eigenvectors_list[q_val][mu_val]

            for atom in range(number_atoms):
                amplitude = phonon_amplitude(temperature, eigenvalue, mass_list[atom])
                initial_phase = random_uniform(0, 2*np.pi)

                vector[atom][0] = vector[atom][0] + (eigenvectors[atom][0] * amplitude * np.cos(initial_phase))
                vector[atom][1] = vector[atom][1] + (eigenvectors[atom][1] * amplitude * np.cos(initial_phase))
                vector[atom][2] = vector[atom][2] + (eigenvectors[atom][2] * amplitude * np.cos(initial_phase))
                
    for atom in range(number_atoms):
        vector[atom][0] = vector[atom][0] / number_q
        vector[atom][1] = vector[atom][1] / number_q
        vector[atom][2] = vector[atom][2] / number_q

    return vector

def get_phonopy(path_file):
    """
    It extracts the number of atoms and q-points, and the eigenvalues and eigenvectors for each q-point,
    from mesh.yaml phonopy output file. It also extracts the atomic masses of the atoms

    Inputs:
        path_file -> path to the mesh.yaml file
    """

    with open(path_file, "r") as file:
        data = yaml.safe_load(file)

    number_phonons = data['natom'] * 3

    number_atoms = data['natom']
    number_q = data['nqpoint']
    eigenvalues = []
    eigenvectors = []
    mass_list = []

    for q_val in range(number_q):
        eigenvalue = []
        for mode in range(number_phonons):
            eigenvalue.append(data['phonon'][q_val]['band'][mode]['frequency'] * 1e12) # to express it in Hz (not in THz)

        eigenvalues.append(eigenvalue)
    
    for q_val in range(number_q):
        eigenvector_mu = []
        for mode in range(number_phonons):
            eigenvector_atom = []
            for atom in range(number_phonons // 3):
                vec_x = data['phonon'][q_val]['band'][mode]['eigenvector'][atom][0][0]
                vec_y = data['phonon'][q_val]['band'][mode]['eigenvector'][atom][1][0]
                vec_z = data['phonon'][q_val]['band'][mode]['eigenvector'][atom][2][0]

                eigenvector_atom.append([vec_x, vec_y, vec_z])

            eigenvector_mu.append(eigenvector_atom)
        
        eigenvectors.append(eigenvector_mu)

    for atom in range(number_phonons // 3):
        mass_list.append(data['points'][atom]['mass'])

    return number_atoms, number_q, eigenvalues, eigenvectors, mass_list

def new_POSCAR(old_path, disp_vec, new_path):
    """
    Applies the find vector displacement to a given POSCAR

    Inputs:
        old_path -> path to the old POSCAR file (IN CARTESIAN COORDINATES!!!!!)
        disp_vec -> generated displacement vector
        new_path -> path to the new distorted structure
    """
    old = open(old_path, 'r')
    new = open(new_path, 'w')

    for _ in range(7):
        line = old.readline()
        new.write(line)

    num_atoms = 0
    for it in range(len(line.split())):
        num_atoms = num_atoms + int(line.split()[it])

    line = old.readline()
    new.write(line)

    for atom in range(num_atoms):
        line = old.readline()

        pos_x = float(line.split()[0]) + disp_vec[atom][0]
        pos_y = float(line.split()[1]) + disp_vec[atom][1]
        pos_z = float(line.split()[2]) + disp_vec[atom][2]

        new.write(f'     {pos_x:.9f}     {pos_y:.9f}     {pos_z:.9f}\n')

    old.close()
    new.close()


### Code
# Read the mesh.yaml file and obtain phonon data
num_atoms, nq, values, vectors, masses = get_phonopy('phonon-data/anharmonic-T200-nac/mesh.yaml')

# Generate a number of Monte Carlo distorted structures and save them in the desired path
path_to_save = 'distorted_structures/'
number_of_structures = 100

for num_struc in range(number_of_structures):
    result = disp_vect(num_atoms, nq, num_atoms*3, values, vectors, masses, 600)
    path_new = path_to_save + 'POSCAR-' + str(num_struc + 1).zfill(4)
    new_POSCAR('POSCAR.vasp', result, path_new)