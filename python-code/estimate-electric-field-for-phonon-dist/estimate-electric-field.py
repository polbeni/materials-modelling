# Pol Benítez Colominas, July 2025
# Universitat Politècnica de Catalunya

# Estimates the electric field necessary to induce a phonon mode distortion with a given amplitude

def effective_born_charges(born_tensor, phonon_eigenvectors):
    """
    Determines the effective Born charges for a given phonon mode

    Inputs:
        born_tensor: Born charges for each atom and mode
        phonon_eigenvectors: phonon eigenvectors for each mode
    """

    effective_born = [0, 0, 0]

    for atom in range(len(phonon_eigenvectors)):
        effective_born[0] = effective_born[0] + ((born_tensor[atom][0][0] * phonon_eigenvectors[atom][0]) + (born_tensor[atom][1][0] * phonon_eigenvectors[atom][1]) + (born_tensor[atom][2][0] * phonon_eigenvectors[atom][2]))
        effective_born[1] = effective_born[1] + ((born_tensor[atom][0][1] * phonon_eigenvectors[atom][0]) + (born_tensor[atom][1][1] * phonon_eigenvectors[atom][1]) + (born_tensor[atom][2][1] * phonon_eigenvectors[atom][2]))
        effective_born[2] = effective_born[2] + ((born_tensor[atom][0][2] * phonon_eigenvectors[atom][0]) + (born_tensor[atom][1][2] * phonon_eigenvectors[atom][1]) + (born_tensor[atom][2][2] * phonon_eigenvectors[atom][2]))
    
    return effective_born

def electric_field(phonon_eigenvalue, displacement, born_charges):
    """
    Estimates the electric field in order to induce a phonon mode with a given amplitude distortion

    Inputs:
        phonon_eigenvalue: energy of the phonon mode
        displacement: displacement of the mode
        born_charges: effective Born charges of the system
    """

    phonon_eigenvalue = phonon_eigenvalue * 0.0041 # in eV

    e_field = [] # in eV/Angstrom

    if born_charges[0] < 1e-6:
        e_field.append(0)
    else:
        e_field.append(((phonon_eigenvalue**2) * displacement[0]) / born_charges[0])
    
    if born_charges[1] < 1e-6:
        e_field.append(0)
    else:
        e_field.append(((phonon_eigenvalue**2) * displacement[1]) / born_charges[1])

    if born_charges[2] < 1e-6:
        e_field.append(0)
    else:
        e_field.append(((phonon_eigenvalue**2) * displacement[2]) / born_charges[2])


    e_field_si = [e_field[0] * 1.602e10, e_field[1] * 1.602e10, e_field[2] * 1.602e10] # in V/m

    return e_field, e_field_si 

##### Ag3SBr #####
print('For Ag3SBr:')
# Values
born_tensor = [[[1.10857995, 0.00000000, 0.00000000], [0.00000000, 1.10857995, 0.00000000], [0.00000000, 0.00000000, 0.93669914]],
               [[0.93669914, 0.00000000, 0.00000000], [0.00000000, 1.10857995, 0.00000000], [0.00000000, 0.00000000, 1.10857995]],
               [[1.10857995, 0.00000000, 0.00000000], [0.00000000, 0.93669914, 0.00000000], [0.00000000, 0.00000000, 1.10857995]],
               [[-1.03052658, 0.00000000, 0.00000000], [0.00000000, -1.03052658, 0.00000000], [0.00000000, 0.00000000, -1.03052658]],
               [[-2.12333247, 0.00000000, 0.00000000], [0.00000000, -2.12333247, 0.00000000], [0.00000000, 0.00000000, -2.12333247]]]

eigenvalue = 0.4016814169 # THz

eigenvectors = [[0.48549225846682, 0, -0.00322143717322],
                [0.48549225846682, 0, 0.00267084351467],
                [-0.58557635448751, 0, 0.00267084351468],
                [-0.37672125956292, 0, -0.00207246050869],
                [-0.20915471808683, 0, -0.00115062498449]]

distortion = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0] # Angstrom


# Determination of the electric field
eff_born = effective_born_charges(born_tensor, eigenvectors)

for dist in distortion:
    _, e_field_si = electric_field(eigenvalue, [dist, 0, 0], eff_born)
    print(f'The electric field necessary to induce a displacement of {dist} Agnstrom is: {e_field_si[0]} V/m')

print('')