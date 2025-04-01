# Pol Benítez Colominas, April 2025
# Universitat Politècnica de Catalunya

# It generates a molecular crystal using PyXtal

from pyxtal import pyxtal
from pyxtal.molecule import pyxtal_molecule

def generate_atom_xyz(atom_type):
    """
    It generates a xyz from a given atom in order to can use the molecular functionality of PyXtal

    Inputs:
        atom_type: atom symbol
    """

    xyz_file = open(f'{atom_type}.xyz', 'w')
    xyz_file.write('1\n')
    xyz_file.write(f'{atom_type} atom\n')
    xyz_file.write(f'{atom_type}   0.0   0.0   0.0')
    xyz_file.close()

#### Inputs ####
dimension = 3
symmetry_group = 1
elements_list = ['Benzene.xyz', 'Adamantane.xyz', 'Na', 'I']
stoichiometry = [1, 1, 2, 2]
################

# Create and array with the different molecules and atoms objects for PyXtal
elements_object = []
for element in elements_list:
    if element.endswith('.xyz'):
        elements_object.append(pyxtal_molecule(element))
    else:
        generate_atom_xyz(element)
        elements_object.append(pyxtal_molecule(element + '.xyz'))

# Generate the phase
xtal = pyxtal(molecular=True)
xtal.from_random(dim=dimension, group=symmetry_group, species=elements_object, numIons=stoichiometry)

print("Successfully generated the molecular crystal cif")

# Save the crystal as a cif
xtal.to_file('molecular_crystal.cif')