# Pol Benítez Colominas, February 2025
# Universitat Politècnica de Catalunya

# Script to run a single shot calculation using ASE calculator and MACE machine-learning interatomic potentials
# The returned values are the system energy, the atoms forces and the system stress tensor

import warnings

from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor

from mace.calculators import mace_mp
import torch
torch.serialization.add_safe_globals([slice])

warnings.simplefilter('ignore')

# Read the system (with VASP structure format)
crystal_structure = Structure.from_file('POSCAR')
ase_adaptor = AseAtomsAdaptor()
atoms = ase_adaptor.get_atoms(crystal_structure)

# Perform the calculation
atoms.calc = mace_mp(model='large', device='cpu')

# Obtain the values
energy = atoms.get_potential_energy()
forces = atoms.get_forces()
stress = atoms.get_stress()

# Save the outputs in files
# Energy (single value)
output_file = open('energy.output', 'w')
output_file.write(f'Energy [eV]: {energy}')
output_file.close()
# Forces (F_x, F_y, F_z, for each atom (in rows))
output_file = open('forces.output', 'w')
output_file.write('# atom      F_x [eV/Ang]      F_y [eV/Ang]      F_z [eV/Ang]\n')
for num_atom in range(len(forces)):
    output_file.write(f'{num_atom + 1}      {forces[num_atom][0]}      {forces[num_atom][1]}      {forces[num_atom][2]}\n')
output_file.close()
# Stress (stress tensor in Voigt order (xx, yy, zz, yz, xz, xy))
output_file = open('stress.output', 'w')
output_file.write('xx [eV/Ang^2]      yy [eV/Ang^2]      zz [eV/Ang^2]      yz [eV/Ang^2]      xz [eV/Ang^2]      xy [eV/Ang^2]\n')
output_file.write(f'{stress[0]}      {stress[1]}      {stress[2]}      {stress[3]}      {stress[4]}      {stress[5]}')
output_file.close()