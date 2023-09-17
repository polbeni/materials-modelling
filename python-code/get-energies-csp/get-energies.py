# Pol Benítez Colominas
# Universitat Politècnica de Catalunya, September 2023

import os
import shutil
import subprocess

num_phases = 15 # number of phases we want to relax
dir_name = 'Ag3SBr' # name of the folder of the material of interest

energy_file_path = dir_name + '/relax-phases/energy_ranking.txt'

phases_relaxed = []

energy_file = open(energy_file_path, "r")

energy_file.readline()
for num in range(num_phases):
    line = energy_file.readline()

    phases_relaxed.append(line.split()[1])  

def get_ground_state(atoms_unit_cell, oszicar_file):
    """
    This functions gets the ground state for a VASP simulation results

    atoms_unit_cell: number of atoms in the unit cell
    oszicar_file: path direction for the OSZICAR file
    """


    file = open(oszicar_file, "r")

    with file as f:
        for line in f:
            pass
        last_line = line

    ground_state = float(last_line.split()[4])

    file.close()

    ground_state_energy_fu = ground_state/float(atoms_unit_cell)

    return ground_state_energy_fu

phases_file_path = dir_name + '/final_phases.txt'
phases_file = open(phases_file_path, "w")

phases_file.write('POSCAR-num   energy per atom M3GNet (eV)   Phase group M3GNet   energy per atom  PBEsol (eV)   Phase group PBEsol \n')

energy_file = open(energy_file_path, "r")

energy_file.readline()

for phase in phases_relaxed:
    poscar_name = phase

    line = energy_file.readline()
    energy_m3gnet = line.split()[2]

    vaspkit_poscar_path = dir_name + '/relax-phases/' + phase + '/vaspkit_POSCAR.log'
    vaspkit_poscar = open(vaspkit_poscar_path, "r")
    for x in range(36):
        vaspkit_poscar.readline()
    line_vaspkit_poscar = vaspkit_poscar.readline()
    num_atoms = line_vaspkit_poscar.split()[2]
    for x in range(11):
        vaspkit_poscar.readline()
    line_vaspkit_poscar = vaspkit_poscar.readline()
    group_m3gnet = line_vaspkit_poscar.split()[1]

    oszicar_path = dir_name + '/relax-phases/' + phase + '/OSZICAR'
    energy_pbesol = get_ground_state(num_atoms, oszicar_path)

    source_file = 'vaspkit-scripts/input_vaspkit_CONTCAR'
    destination_file = dir_name + '/relax-phases/' + phase
    shutil.copy(source_file, destination_file)
    
    source_file = 'vaspkit-scripts/vaspkit-CONTCAR.sh'
    shutil.copy(source_file, destination_file)

    bash_path = dir_name + '/relax-phases/' + phase + '/vaspkit-CONTCAR.sh'
    subprocess.run(["bash", bash_path])

    vaspkit_contcar_path = dir_name + '/relax-phases/' + phase + '/vaspkit_CONTCAR.log'
    vaspkit_contcar = open(vaspkit_contcar_path, "r")
    for x in range(48):
        vaspkit_contcar.readline()
    line_vaspkit_contcar = vaspkit_contcar.readline()
    group_pbseol = line_vaspkit_contcar.split()[1]

    phases_file.write(f'{poscar_name}   {energy_m3gnet}          {group_m3gnet}      {energy_pbesol}               {group_pbseol} \n')

    vaspkit_poscar.close()
    vaspkit_contcar.close()

energy_file.close()
phases_file.close()