from pymatgen.core import Structure

import warnings

from m3gnet.models import Relaxer
from pymatgen.core import Lattice, Structure

for category in (UserWarning, DeprecationWarning):
    warnings.filterwarnings("ignore", category=category, module="tensorflow")


# Import the POSCAR file
crystal_structure = Structure.from_file("POSCAR")

# Create the ionic relaxation method
relaxer = Relaxer()

relax_results = relaxer.relax(crystal_structure, verbose=True)

# Show the results
final_structure = relax_results['final_structure']
final_energy_per_atom = float(relax_results['trajectory'].energies[-1] / len(crystal_structure))

print("The relaxed structure is: ", final_structure)
print(f"Final energy is {final_energy_per_atom:.3f} eV/atom")

# Save the structure in CONTCAR type file
final_structure.to(filename='CONTCAR', fmt='poscar')