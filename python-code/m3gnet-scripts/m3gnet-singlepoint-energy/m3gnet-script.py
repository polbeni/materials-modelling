from pymatgen.core import Structure

import warnings

from pymatgen.io.ase import AseAtomsAdaptor

import matgl
from matgl.ext.ase import PESCalculator


warnings.simplefilter("ignore") # To suppress warnings for clearer output


# Import the POSCAR file
crystal_structure = Structure.from_file("POSCAR")

# Load the M3GNet model
pot = matgl.load_model("M3GNet-MP-2021.2.8-PES")

# Create the ionic relaxation method
calc = PESCalculator(potential=pot)

# Create ase atoms object
ase_adaptor = AseAtomsAdaptor()
atoms = ase_adaptor.get_atoms(crystal_structure)

# Compute the energy and print the results
atoms.set_calculator(calc)
print(f"The calculated potential energy is {atoms.get_potential_energy():.3f} eV.")