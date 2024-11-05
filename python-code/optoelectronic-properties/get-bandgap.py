# Pol Benítez Colominas, November 2024
# University of Cambridge and Universitat Politècnica de Catalunya

# Determine the band gap from vasprun.xml file using pymatgen module

from pymatgen.io.vasp.outputs import Vasprun

# Load the vasprun.xml file
vasprun = Vasprun("vasprun.xml")

# Extract the band structure
band_structure = vasprun.get_band_structure()

# Calculate the band gap
band_gap = band_structure.get_band_gap()

# Print the results
print(f"Band gap type: {band_gap['direct']}")
print(f"Band gap energy (eV): {band_gap['energy']}")