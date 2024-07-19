# Pol Benítez Colominas, July 2024
# Universitat Politècnica de Catalunya

# Compute the surface area of a slab

from pymatgen.core import Structure
import numpy as np

def get_area(path_struc, lattice_param_1, lattice_param_2):
    """
    Compute the area of a unit cell surface

    Inputs:
        path_struc: path to POSCAR file
        lattice_param_1, lattice_param_2: lattice parameters that create the surface (0->x, 1->y, 2->z)
    """

    structure = Structure.from_file(path_struc)

    lattice = structure.lattice
    a = lattice.matrix[lattice_param_1]
    b = lattice.matrix[lattice_param_2]

    cross_product = np.cross(a, b)

    area = np.linalg.norm(cross_product)

    return area

print(f"The surface area perpendicular to the z-axis for (100) slab is: {get_area('Ag3SBr/100.vasp', 0, 1)} Å^2")
print(f"The surface area perpendicular to the z-axis for (110) slab is: {get_area('Ag3SBr/110.vasp', 0, 1)} Å^2")
print(f"The surface area perpendicular to the z-axis for (111) slab is: {get_area('Ag3SBr/111.vasp', 0, 1)} Å^2")