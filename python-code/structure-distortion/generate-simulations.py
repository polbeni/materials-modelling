# Pol Benítez Colominas, September 2024
# Universitat Politècnica de Catalunya and University of Cambridge

# This code generates the distorted POSCAR files for the desired eigenvectors at gamma-point in
# order to study LO-TO splitting effect in e-p band gap

import shutil
import yaml
import subprocess

def execute_claudi(path):
    """
    Executes Claudi program to distort the structure

    Inputs:
        path: path where the calculation is performed
    """
    #exe_path = './displ.exe'

    input1 = '3' # number of different species in the unit cell
    input2 = '5' # total number of atoms in the unit cell
    input3 = '0.4' # amplitude displacment in angstroms
    input4 = '0' # to include normalization of eigenvectors or not (0: No, 1: Yes)
    inputs = f'{input1}\n{input2}\n{input3}\n{input4}\n'

    result = subprocess.run([path], input=inputs, text=True, capture_output=True)

    return

# cartesian positions of the atoms
cartesian_positions = [
    "  2.396690  2.396690  0.000000   ",
    "  2.396690  0.000000  2.396690   ",
    "  0.000000  2.396690  2.396690   ",
    "  2.396690  2.396690  2.396690   ",
    "  0.000000  0.000000  0.000000   "
]

def get_eigenvectors(file_path, band_number, atom_number):
    """
    Returns the eigenvectors for a given band and given atom, at the gamma-point, using phonopy band.yaml output file

    Inputs:
        file_path: path to the band.yaml file
        band_number: number of the phonon band
        atom_number: number of the atom
    """

    with open(file_path, "r") as file:
        data = yaml.safe_load(file)
    phonon = data['phonon'][0]['band'][band_number]['eigenvector'][atom_number]
    
    eigenvector_string = f'{str(phonon[0][0])}  {str(phonon[1][0])}  {str(phonon[2][0])}'

    return eigenvector_string

def generate_VECTORS(file_path, band_number, output_path):
    """
    Generates the VECTORS file for a given band

    Inputs:
        file_path: path to the band.yaml file
        band_number: number of the phonon band
        output_path: path to save the output
    """
    VECTORS = open(output_path, 'w')
    
    for x in range(5): # iterate over all the atoms in the unit cell
        results = get_eigenvectors(file_path, band_number, x)
        VECTORS.write(f'{cartesian_positions[x]}{results}\n')

    VECTORS.close()

phonons_calc = ['phonon-01', 'phonon-02', 'phonon-03', 'phonon-04', 'phonon-05', 
                'phonon-06', 'phonon-07', 'phonon-08', 'phonon-09', 'phonon-10', 
                'phonon-11', 'phonon-12', 'phonon-13', 'phonon-14', 'phonon-15']

# no-nac
bands_file = '../../harmonic-2x2x2/no-nan/eigenvectors/band.yaml'
dir_name = 'no-nac/'

num_band = 0
for band in phonons_calc:
    print(band)

    generate_VECTORS(bands_file, num_band, 'VECTORS')
    execute_claudi('./displ.exe')

    # POSCAR
    source = 'POSCARnew'
    destination = dir_name + band + '/POSCAR'
    shutil.copy(source, destination)

    # KPOINTS
    source = 'save/KPOINTS'
    destination = dir_name + band + '/KPOINTS'
    shutil.copy(source, destination)

    # INCAR
    source = 'save/INCAR'
    destination = dir_name + band + '/INCAR'
    shutil.copy(source, destination)

    # POTCAR
    source = 'save/POTCAR'
    destination = dir_name + band + '/POTCAR'
    shutil.copy(source, destination)

    # run.sh
    source = 'save/run.sh'
    destination = dir_name + band + '/run.sh'
    shutil.copy(source, destination)

    num_band = num_band + 1

# nac
bands_file = '../../harmonic-2x2x2/nan-wang/eigenvectors/band.yaml'
dir_name = 'nac/'

num_band = 0
for band in phonons_calc:
    print(band)
    
    generate_VECTORS(bands_file, num_band, 'VECTORS')
    execute_claudi('./displ.exe')

    # POSCAR
    source = 'POSCARnew'
    destination = dir_name + band + '/POSCAR'
    shutil.copy(source, destination)

    # KPOINTS
    source = 'save/KPOINTS'
    destination = dir_name + band + '/KPOINTS'
    shutil.copy(source, destination)

    # INCAR
    source = 'save/INCAR'
    destination = dir_name + band + '/INCAR'
    shutil.copy(source, destination)

    # POTCAR
    source = 'save/POTCAR'
    destination = dir_name + band + '/POTCAR'
    shutil.copy(source, destination)

    # run.sh
    source = 'save/run.sh'
    destination = dir_name + band + '/run.sh'
    shutil.copy(source, destination)

    num_band = num_band + 1