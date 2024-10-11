# Pol Benítez Colominas, October 2024
# University of Cambridge and Universitat Politècnica de Catalunya

# Computes the dielectric constant for an anisotropic dielectric tensor following https://doi.org/10.1038/s41524-023-01083-8
# and computes the mean value for the thermal structures

import numpy as np

def get_dielectric(OUTCAR_file):
    """
    Reads the OUTCAR file from VASP simulation and returns the different dielectric tensors

    Inputs:
        OUTCAR_file -> path to the OUTCAR file
    Outputs:
        dielectric_el -> returns the electronic part of the dielectric tensor
                         in VASP: MACROSCOPIC STATIC DIELECTRIC TENSOR
        dielectric_ion -> returns the ionic part of the dielectric tensor
                          in VASP: MACROSCOPIC STATIC DIELECTRIC TENSOR IONIC CONTRIBUTION
    """

    with open(OUTCAR_file, 'r') as file:
        searc_pattern = 'MACROSCOPIC STATIC DIELECTRIC TENSOR (including local field effects in DFT)'

        num_lines = 1
        count = True
        for line in file:
            if searc_pattern in line:
                count = False
            
            if count == True:
                num_lines = num_lines + 1

    OUTCAR = open(OUTCAR_file, 'r')
    for x in range(num_lines):
        OUTCAR.readline()
    line = OUTCAR.readline()
    line = OUTCAR.readline()
    el11 = float(line.split()[0])
    el12 = float(line.split()[1])
    el13 = float(line.split()[2])
    line = OUTCAR.readline()
    el21 = float(line.split()[0])
    el22 = float(line.split()[1])
    el23 = float(line.split()[2])
    line = OUTCAR.readline()
    el31 = float(line.split()[0])
    el32 = float(line.split()[1])
    el33 = float(line.split()[2])

    dielectric_el = np.matrix([[el11, el12, el13],
                               [el21, el22, el23],
                               [el31, el32, el33]])
    OUTCAR.close()

    with open(OUTCAR_file, 'r') as file:
        searc_pattern = 'MACROSCOPIC STATIC DIELECTRIC TENSOR IONIC CONTRIBUTION'

        num_lines = 1
        count = True
        for line in file:
            if searc_pattern in line:
                count = False
            
            if count == True:
                num_lines = num_lines + 1

    OUTCAR = open(OUTCAR_file, 'r')
    for x in range(num_lines):
        OUTCAR.readline()
    line = OUTCAR.readline()
    line = OUTCAR.readline()
    el11 = float(line.split()[0])
    el12 = float(line.split()[1])
    el13 = float(line.split()[2])
    line = OUTCAR.readline()
    el21 = float(line.split()[0])
    el22 = float(line.split()[1])
    el23 = float(line.split()[2])
    line = OUTCAR.readline()
    el31 = float(line.split()[0])
    el32 = float(line.split()[1])
    el33 = float(line.split()[2])

    dielectric_ion = np.matrix([[el11, el12, el13],
                                [el21, el22, el23],
                                [el31, el32, el33]])
    OUTCAR.close()

    return dielectric_el, dielectric_ion

def dielectric_ct(dielectric_tensor):
    """
    It returns the dielectric constant for the anisotropic dielectric tensor

    Inputs:
        dielectric_tensor -> dielectric tensor we want to get the effective value
    """

    inverse_matrix = np.linalg.inv(dielectric_tensor)
    trace = np.trace(inverse_matrix) / 3
    effective = 1 / trace

    return effective

diel1, diel2 = get_dielectric('results/sim-01/OUTCAR')
print(diel1, dielectric_ct(diel1))
print(diel2, dielectric_ct(diel2))