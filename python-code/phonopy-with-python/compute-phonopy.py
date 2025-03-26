# Pol Benítez Colominas, March 2025
# Universitat Politècnica de Catalunya

# Use phonopy to compute phonons and thermal properties
# The forces are determined using MLIP models (such as MACE)

import numpy as np

from pymatgen.core import Structure, Lattice
from pymatgen.io.ase import AseAtomsAdaptor

from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms

from mace.calculators import mace_mp

import warnings
warnings.filterwarnings("ignore")

# read the structure with pymatgen (POSCAR or cif file)
structure = Structure.from_file('POSCAR')

# generate the unit cell in phonopy format
unit_cell = PhonopyAtoms(symbols=[site.species_string for site in structure],
                         scaled_positions=structure.frac_coords,
                         cell=structure.lattice.matrix)

# generate supercell and initialize phonopy
supercell_matrix = [[3, 0, 0], [0, 3, 0], [0, 0, 3]]

phonon = Phonopy(unit_cell, supercell_matrix)

# generate displacementes
distance = 0.01 # this is the default value in phonopy
phonon.generate_displacements(distance=distance)
supercells = phonon.supercells_with_displacements

print(f"Supercell contains {phonon.supercell.get_number_of_atoms()} atoms")
print(f"Generated {len(phonon.supercells_with_displacements)} displacements")

# transform phonopy distorted supercell strucutres to pymatgen structure
def phonopy_to_pymatgen(phonopy_atoms):
    """
    Transforms phonopy structure into pymatgen structure

    Inputs:
        phonopy_atoms: object with phonopy structure
    """

    lattice = Lattice(phonopy_atoms.cell)
    species = phonopy_atoms.symbols
    coords = phonopy_atoms.scaled_positions
    
    return Structure(lattice=lattice, species=species, coords=coords)

# compute the forces for each distorted supercell using MACE
forces = []

it_supercell = 1
for supercell in supercells:
    print(f'Computing forces for supercell {it_supercell} of a total of {len(supercells)} supercells')
    supercell_structure = phonopy_to_pymatgen(supercell)

    ase_adaptor = AseAtomsAdaptor()
    atoms = ase_adaptor.get_atoms(supercell_structure)
    atoms.calc = mace_mp(model='large', device='cpu')

    forces_struc = atoms.get_forces()
    
    forces.append(forces_struc)

    print('')
    it_supercell = it_supercell + 1
    
# generate the force constants
phonon.forces = forces
phonon.produce_force_constants()

# calculate phonon properties (DOS)
mesh = [30, 30, 30]
phonon.set_mesh(mesh, is_mesh_symmetry=True)
phonon.set_total_DOS()

np.savetxt('phonon_dos.dat', np.array(phonon.get_total_DOS()).T, header='Frequency(THz) DOS')

# count the proportion of imaginary freq
def proportion_imaginary(dos):
    """
    Returns the proportion of imaginary phonon modes in the phonon DOS

    Inputs:
        dos: phonon DOS (density of states)
    """
    prop_imag = 0
    prop_real = 0

    for state in range(len(dos[0])):
        if dos[0][state] < 0:
            prop_imag = prop_imag + dos[1][state]
        else:
            prop_real = prop_real + dos[1][state]
    
    return prop_imag / (prop_imag + prop_real)

prop_imaginary_phonon = proportion_imaginary(np.array(phonon.get_total_DOS()))

print(f'The proportion of imaginary phonon modes with respect real modes is: {prop_imaginary_phonon}')

# compute thermal properties
phonon.set_thermal_properties(t_max=1000, t_min=0, t_step=10)
temperatures, free_energy = phonon.get_thermal_properties()[:2]

np.savetxt('thermal_properties.dat', np.column_stack((temperatures, free_energy)), header='Temperature(K) Free_energy(kJ/mol)')