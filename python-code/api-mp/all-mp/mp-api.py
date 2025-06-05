# Pol Benítez Colominas, June 2025
# Universitat Politècnica de Catalunya

# Retrieve all the materials in the Materials Project database

import os
import shutil

from mp_api.client import MPRester
from pymatgen.io.vasp.inputs import Poscar

# Your Materials Project API key
api_key = "your_api_key"

with MPRester(api_key) as mpr:
    docs = mpr.materials.summary.search(
        fields=["material_id", "structure", "formula_pretty", "elements"]
    )

# Output format for structures
file_fmt = 'cif'

# Create a directory for saving structure files
output_dir = 'all_structures'
if os.path.exists(output_dir):
    shutil.rmtree(output_dir)
os.mkdir(output_dir)

# Open a file to record the information
info_file = open(f'{output_dir}/materials_info.txt', 'w')
info_file.write('Material ID     Formula     Number of Species\n')

# Process each material
for material in docs:
    material_id = material.material_id
    formula = material.formula_pretty
    elements = material.elements
    num_species = len(elements)

    # Save info
    info_file.write(f'{material_id}     {formula}     {num_species}\n')

    # Save structure
    structure = material.structure
    if file_fmt.lower() == 'poscar':
        poscar = Poscar(structure)
        file_name = f'{output_dir}/{material_id}.poscar'
        poscar.write_file(file_name)
    elif file_fmt.lower() == 'cif':
        file_name = f'{output_dir}/{material_id}.cif'
        structure.to(fmt="cif", filename=file_name)
    else:
        print(f'Unsupported file format: {file_fmt}')

info_file.close()
