# Pol Benítez Colominas, February 2025
# Universitat Politècnica de Catalunya

# Apply harmonic phonon dispersion at given T
# It just consider the phonon modes at gamma

import yaml
import numpy as np

### Physical constants
reduced_planck_ct = 6.5821e-16 # eV·s
boltzmann_ct = 8.6173e-5 # eV·K^-1

### Functions
def bose_einstein(temperature, phonon_energy):
    """
    Computes the Bose-Einstein distribution of a phonon population with a given energy

    Inputs:
        temperature -> temperature (in K)
        phonon_energy -> energy of the given phonon (in s^-1 (Hz))
    """

    distribution = 1 / (np.exp((reduced_planck_ct * 2 * np.pi * phonon_energy) / (boltzmann_ct * temperature)) - 1) 

    return distribution

def phonon_amplitude(temperature, phonon_energy, atom_mass):
    """
    Computes the phonon amplitude for a given temperature and phonon mode

    Inputs:
        temperature -> temperature (in K)
        phonon_energy -> energy of the given phonon (in s^-1 (Hz))
        atom_mass -> atomic mass of the given atom
    """

    dist = bose_einstein(temperature, phonon_energy)

    h_si = reduced_planck_ct * 6.24e-18 # in J·s
    mass = atom_mass * 1.660539e-27 # in Kg

    amplitude = np.sqrt((h_si / (2 * mass * 2 * np.pi * phonon_energy)) * (1 + 2*dist))

    amplitude = amplitude * 1e10 # in angstrom

    return amplitude

def get_phonopy(path_file):
    """
    It extracts the number of atoms and the eigenvalues and eigenvectors at gamma-point,
    from band.yaml phonopy output file. It also extracts the atomic masses of the atoms
    IMPORTANT: it assumes that gamma point is the first q-pont in the band.yaml file
    change the k-path to start at gamma otherwise

    Inputs:
        path_file -> path to the band.yaml file
    """

    with open(path_file, "r") as file:
        data = yaml.safe_load(file)

    number_phonons = data['natom'] * 3

    number_atoms = data['natom']
    number_q = 1
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

    return number_atoms, eigenvalues, eigenvectors, mass_list

def compute_norm(eigenvectors, num_atoms, acoustic):
    """
    It computes the norm of a given set of eigenvectors for all the phonon modes
    neglecting the acoustics since we are at gamma

    Inputs:
        eigenvectors -> phonon eigenvectors at the gamma-point
        num_atoms -> number of atoms in the structure
        acoustic -> boolean indicating if the list of eigenvectors contain acoustics (True) or not (False)
    """

    phonon_modes = (num_atoms * 3) - 3

    norm = 0

    for mode in range(phonon_modes):
        for atom in range(num_atoms):
            for coord in range(3):
                if acoustic == True:
                    norm = norm + ((eigenvectors[mode + 3][atom][coord])**2)
                elif acoustic == False:
                    norm = norm + ((eigenvectors[mode][atom][coord])**2)

    norm = np.sqrt(norm)

    return norm

def normalize_eigenvectors(eigenvectors, num_atoms):
    """
    Computes the normalization constant and returns renormalized eigenvectors

    Inputs:
        eigenvectors -> phonon eigenvectors of a given phonon mode at the gamma-point
        num_atoms -> number of atoms in the structure
    """

    normalization = 1 / compute_norm(eigenvectors, num_atoms, True)

    phonon_modes = (num_atoms * 3) - 3

    new_vectors = []

    for mode in range(phonon_modes):
        mode_vectors = []
        for atom in range(num_atoms):
            atom_coords = [(eigenvectors[mode + 3][atom][0]) * normalization,
                           (eigenvectors[mode + 3][atom][1]) * normalization,
                           (eigenvectors[mode + 3][atom][2]) * normalization]
                
            mode_vectors.append(atom_coords)
        
        new_vectors.append(mode_vectors)

    return new_vectors

def disp_mode(eigenvectors, num_atoms, temperature, phonon_energy, masses_list):
    """
    It generates the displacement vector for a given phonon mode and a given temperature
    It just considers the eigenvectors at gamma-point

    Inputs:
        eigenvectors -> phonon eigenvectors for the desired mode (normalized)
        num_atoms -> number of atoms in the structure
        temperature -> temperature (in K)
        phonon_energy -> energy of the given phonon (in s^-1 (Hz))
        masses_list -> list with the atomic masses of the atoms
    """

    disp_vect = []

    for atom in range(num_atoms):
        amplitude = phonon_amplitude(temperature, phonon_energy, masses_list[atom])

        atom_vect = [eigenvectors[atom][0] * amplitude,
                     eigenvectors[atom][1] * amplitude,
                     eigenvectors[atom][2] * amplitude]
        
        disp_vect.append(atom_vect)

    return disp_vect


### Code
# Get the data at gamma-point
num_atoms, values, vectors, masses = get_phonopy('band.yaml')

# Renormalize the eigenvectors
norm_vectors = normalize_eigenvectors(vectors[0], num_atoms)

# Get the final displacement vectors
disp_final = []
for x in range(12):
        disp_final.append(disp_mode(norm_vectors[x], num_atoms, 1000, values[0][x + 3], masses))

print(disp_final)