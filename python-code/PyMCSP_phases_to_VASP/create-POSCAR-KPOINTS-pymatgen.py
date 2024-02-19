# Pol Benítez Colominas
# Universitat Politècnica de Catalunya, September 2023

import os
import shutil

from pymatgen.io.vasp import Poscar

num_phases = 20 # number of phases we want to relax
material = 'Ag6SSeI2'
dir_name = '../' + material # name of the folder of the material of interest

energy_file_path = dir_name + '/structure_files/initial_structures/relaxed_structures/energy_ranking.txt'

phases_to_relax = []

energy_file = open(energy_file_path, "r")

energy_file.readline()
for num in range(num_phases):
    line = energy_file.readline()

    phases_to_relax.append(line.split()[1])


if os.path.exists('relax-phases/' + material):
    shutil.rmtree('relax-phases/' + material)
os.mkdir('relax-phases/' + material)

source_file = dir_name + '/structure_files/initial_structures/relaxed_structures/energy_ranking.txt' 
destination_file = 'relax-phases/' + material
shutil.copy(source_file, destination_file)

for phase in phases_to_relax:
    phase_new_path = 'relax-phases/'  + material + '/' + phase
    os.mkdir(phase_new_path)

    source_file = dir_name + '/structure_files/initial_structures/relaxed_structures/' + phase
    destination_file = phase_new_path + '/POSCAR'
    shutil.copy(source_file, destination_file)

energy_file.close()

def number_KPOINTS(lattice_ct_a, lattice_ct_b, lattice_ct_c):
    """
    Gives the number of KPOINTS for each direction

    Inputs:
        lattice_ct_a, lattice_ct_b, lattice_ct_c: lattice constants for the three directions
    """

    lattice_fu = 5    # size im angstroms
    kpoints_fu = 6    # kpoints for lattice_fu 

    kpoints_a = int(kpoints_fu*(lattice_fu/lattice_ct_a))
    kpoints_b = int(kpoints_fu*(lattice_fu/lattice_ct_b))
    kpoints_c = int(kpoints_fu*(lattice_fu/lattice_ct_c))

    return kpoints_a, kpoints_b, kpoints_c

for phase in phases_to_relax:
    phase_new_path = 'relax-phases/' + material + '/' + phase
    KPOINTS_path = phase_new_path + '/KPOINTS'

    KPOINTS = open(KPOINTS_path, "w")
    KPOINTS.write('Automatic mesh \n')
    KPOINTS.write('0 \n')
    KPOINTS.write('G \n')

    POSCAR_path = phase_new_path + '/POSCAR'
    poscar = Poscar.from_file(POSCAR_path)

    lattice_parameters = poscar.structure.lattice.abc

    
    latt_a = lattice_parameters[0]
    latt_b = lattice_parameters[1]
    latt_c = lattice_parameters[2]

    kpoinpt_a, kpoint_b, kpoint_c = number_KPOINTS(latt_a, latt_b, latt_c)

    KPOINTS.write(f' {kpoinpt_a}   {kpoint_b}   {kpoint_c} \n')

    KPOINTS.write(' 0.  0.  0.')
    KPOINTS.close()

path_scripts = 'relax-phases/' + material +'/'

cp_script = open(path_scripts + 'cp_script.sh', "w")
cp_script.write('#!/bin/bash \n')
cp_script.write(' \n')
phases_line = 'phases=('
for phase in phases_to_relax:
    if phase == phases_to_relax[-1]:
        phases_line = phases_line + '"' + phase + '"'
    else:
        phases_line = phases_line + '"' + phase + '" '
phases_line = phases_line + ') \n'
cp_script.write(phases_line)
cp_script.write(' \n')
cp_script.write('for dir in "${phases[@]}"; do \n')
cp_script.write('cp save/* $dir \n')
cp_script.write('done')
cp_script.close()
os.system("chmod +x " + path_scripts + "cp_script.sh")

send_jobs = open(path_scripts + 'send_jobs.sh', "w")
send_jobs.write('#!/bin/bash \n')
send_jobs.write(' \n')
phases_line = 'phases=('
for phase in phases_to_relax:
    if phase == phases_to_relax[-1]:
        phases_line = phases_line + '"' + phase + '"'
    else:
        phases_line = phases_line + '"' + phase + '" '
phases_line = phases_line + ') \n'
send_jobs.write(phases_line)
send_jobs.write(' \n')
send_jobs.write('for dir in "${phases[@]}"; do \n')
send_jobs.write('cd $dir \n')
send_jobs.write('sbatch run.sh \n')
send_jobs.write('cd .. \n')
send_jobs.write('done')
send_jobs.close()
os.system("chmod +x " + path_scripts + "send_jobs.sh")

check_simulations = open(path_scripts + 'check_simulations.sh', "w")
check_simulations.write('#!/bin/bash \n')
check_simulations.write(' \n')
phases_line = 'phases=('
for phase in phases_to_relax:
    if phase == phases_to_relax[-1]:
        phases_line = phases_line + '"' + phase + '"'
    else:
        phases_line = phases_line + '"' + phase + '" '
phases_line = phases_line + ') \n'
check_simulations.write(phases_line)
check_simulations.write(' \n')
check_simulations.write('for dir in "${phases[@]}"; do \n')
check_simulations.write('echo simulation $dir \n')
check_simulations.write('tail -n 1 $dir/OUTCAR \n')
check_simulations.write('done')
check_simulations.close()
os.system("chmod +x " + path_scripts + "check_simulations.sh")