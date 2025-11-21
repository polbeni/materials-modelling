# Pol Benítez Colominas, November 2025
# The University of Tokyo and Universitat Politècnica de Catalunya

# Uses FC file from phonopy calculations to generate temperature displacements

import os
import numpy as np

from phonopy import Phonopy
from phonopy.interface.calculator import read_crystal_structure
from phonopy.interface.vasp import write_vasp
from phonopy.file_IO import parse_FORCE_CONSTANTS, parse_FORCE_SETS


def generate_random_displacements(
    poscar_file='POSCAR',
    force_file='FORCE_CONSTANTS',
    force_file_type='FORCE_SETS',
    supercell_matrix=[[2, 0, 0], [0, 2, 0], [0, 0, 2]],
    temperature=100,
    n_snapshots=100,
    fc_symmetry=True,
    output_dir='displaced_structures',
    primitive_matrix='auto'
):
    """
    Generate random thermal displacements using Phonopy API
    
    Parameters:
    -----------
    poscar_file : str
        Path to POSCAR file (unit cell)
    force_file : str
        Path to FORCE_CONSTANTS file
    force_file_type : str
        Type of force file: 'FORCE_CONSTANTS' or 'FORCE_SETS'
    supercell_matrix : list or array
        Supercell matrix (DIM parameter)
    temperature : float
        Temperature for random displacements in Kelvin
    n_snapshots : int
        Number of displaced structures to generate
    fc_symmetry : bool
        Apply symmetry to force constants
    output_dir : str
        Directory to save displaced structures
    primitive_matrix : str or array
        Primitive cell matrix ('auto' or explicit matrix)
    """
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Read unit cell structure
    print(f"Reading structure from {poscar_file}")
    unitcell, _ = read_crystal_structure(poscar_file, interface_mode='vasp')
    
    # Create Phonopy object
    print(f"Creating Phonopy object with supercell matrix: {supercell_matrix}")
    phonon = Phonopy(
        unitcell,
        supercell_matrix=supercell_matrix,
        primitive_matrix=primitive_matrix
    )
    
    # Load forces depending on file type
    if force_file_type.upper() == 'FORCE_CONSTANTS':
        print(f"Reading force constants from {force_file}")
        force_constants = parse_FORCE_CONSTANTS(filename=force_file)
        phonon.set_force_constants(force_constants)
        
    elif force_file_type.upper() == 'FORCE_SETS':
        print(f"Reading force sets from {force_file}")
        force_sets = parse_FORCE_SETS(filename=force_file)
        phonon.dataset = force_sets
        
        # Produce force constants from force sets
        print("Producing force constants from force sets")
        phonon.produce_force_constants()
        
    else:
        raise ValueError(f"Unknown force_file_type: {force_file_type}. Use 'FORCE_CONSTANTS' or 'FORCE_SETS'")
    
    # Apply symmetry if requested
    if fc_symmetry:
        print("Applying symmetry to force constants")
        phonon.symmetrize_force_constants()
    
    # Generate random displacements
    print(f"\nGenerating {n_snapshots} random displacements at T = {temperature} K")
    phonon.generate_displacements(
        number_of_snapshots=n_snapshots,
        random_seed=None,  # Use None for random seed, or set an integer for reproducibility
        temperature=temperature,
        cutoff_frequency=None  # Can set a cutoff to ignore low-frequency modes
    )
    
    # Get the supercells with displacements
    supercells = phonon.supercells_with_displacements
    
    if supercells is None or len(supercells) == 0:
        print("Error: No displacements generated!")
        return
    
    print(f"Generated {len(supercells)} displaced supercells")
    
    # Save each displaced structure
    for i, supercell in enumerate(supercells):
        output_file = os.path.join(output_dir, f'POSCAR-{i+1:03d}')
        write_vasp(output_file, supercell)
        if (i + 1) % 10 == 0:
            print(f"  Saved {i+1}/{len(supercells)} structures")
    
    print(f"\nAll displaced structures saved to {output_dir}/")
    print(f"Files: POSCAR-001 to POSCAR-{len(supercells):03d}")
    
    # Save a summary file
    summary_file = os.path.join(output_dir, 'displacement_info.txt')
    with open(summary_file, 'w') as f:
        f.write(f"Random Displacement Generation Summary\n")
        f.write(f"=" * 50 + "\n\n")
        f.write(f"Input POSCAR: {poscar_file}\n")
        f.write(f"Force file: {force_file} ({force_file_type})\n")
        f.write(f"Supercell matrix: {supercell_matrix}\n")
        f.write(f"Primitive matrix: {primitive_matrix}\n")
        f.write(f"Temperature: {temperature} K\n")
        f.write(f"Number of snapshots: {n_snapshots}\n")
        f.write(f"FC symmetry applied: {fc_symmetry}\n")
        f.write(f"Number of atoms in supercell: {len(supercell.positions)}\n")
    
    print(f"Summary saved to {summary_file}")
    
    return phonon


temp_list = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600]
output_dir = 'output_struc'

for temp in temp_list:
    phonon = generate_random_displacements(poscar_file='POSCAR',
                                        force_file='FORCE_CONSTANTS',
                                        force_file_type='FORCE_CONSTANTS',
                                        supercell_matrix=[[4, 0, 0], [0, 4, 0], [0, 0, 4]],
                                        temperature=temp,
                                        n_snapshots=10,
                                        fc_symmetry=True,
                                        output_dir=f'{output_dir}/{str(temp)}',
                                        primitive_matrix='auto')