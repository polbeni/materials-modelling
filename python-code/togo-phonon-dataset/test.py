from pymatgen.io.vasp import Poscar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

structure = Poscar.from_file('POSCAR1').structure

def check_materials(strcuture, name_structure):
    """
    Check if a material from Togo phonon data set is interesting or not
    Here we are interested in three conditions:
        - Centrosymmetry
        - Large imaginary polar phonon
        - Band gap between 1.5 and 4 eV (wide-bandgap semiconductor)

    Inputs:
        structure: structure from Togo dataset to study
        name_structure: name of the structure in Togo dataset
    """

    

    if (centrosymmetric == True) and (imaginary_phonon == True) and (desired_bandgap == True):

spacegroup_analyzer = SpacegroupAnalyzer(structure)

# Check if the space group contains inversion symmetry
is_centrosymmetric = spacegroup_analyzer.is_laue()

print(spacegroup_analyzer.is_laue())

if is_centrosymmetric:
    print("The structure is centrosymmetric.")
else:
    print("The structure is not centrosymmetric.")