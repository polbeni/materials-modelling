# Pol Benítez Colominas
# Universitat Politècnica de Catalunya, September 2023

import os
import shutil
import subprocess

num_phases = 15 # number of phases we want to relax
dir_name = 'Ag3SBr' # name of the folder of the material of interest

energy_file_path = dir_name + '/structure_files/initial_structures/relaxed_structures/energy_ranking.txt'

phases_to_relax = []

energy_file = open(energy_file_path, "r")

energy_file.readline()
for num in range(num_phases):
    line = energy_file.readline()

    phases_to_relax.append(line.split()[1])

if os.path.exists('relax-phases'):
    shutil.rmtree('relax-phases')
os.mkdir('relax-phases')

for phase in phases_to_relax:
    phase_new_path = 'relax-phases/' + phase
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

    kpoints_a = int(8*(5/lattice_ct_a))
    kpoints_b = int(8*(5/lattice_ct_b))
    kpoints_c = int(8*(5/lattice_ct_c))

    return kpoints_a, kpoints_b, kpoints_c

for phase in phases_to_relax:
    phase_new_path = 'relax-phases/' + phase
    KPOINTS_path = phase_new_path + '/KPOINTS'

    KPOINTS = open(KPOINTS_path, "w")
    KPOINTS.write('Automatic mesh \n')
    KPOINTS.write('0 \n')
    KPOINTS.write('G \n')

    source_file = 'vaspkit-scripts/input_vaspkit_POSCAR'
    destination_file = phase_new_path 
    shutil.copy(source_file, destination_file)
    
    source_file = 'vaspkit-scripts/vaspkit-POSCAR.sh'
    shutil.copy(source_file, destination_file)

    bash_path = phase_new_path + '/vaspkit-POSCAR.sh'
    subprocess.run(["bash", bash_path])

    vaspkit_log_path = phase_new_path + '/vaspkit_POSCAR.log'
    vaspkit_log = open(vaspkit_log_path, "r")
    for x in range(42):
        vaspkit_log.readline()
    line_log = vaspkit_log.readline()
    latt_a = float(line_log.split()[2])
    latt_b = float(line_log.split()[3])
    latt_c = float(line_log.split()[4])
    vaspkit_log.close()

    kpoinpt_a, kpoint_b, kpoint_c = number_KPOINTS(latt_a, latt_b, latt_c)

    KPOINTS.write(f' {kpoinpt_a}   {kpoint_b}   {kpoint_c} \n')

    KPOINTS.write(' 0.  0.  0.')
    KPOINTS.close()

cp_script = open('cp_script.sh', "w")
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
cp_script.write('cp save/* relax-phases/$dir \n')
cp_script.write('done')
cp_script.close()
os.system("chmod +x cp_script.sh")

send_jobs = open('send_jobs.sh', "w")
send_jobs.write('#!/bin/bash \n')
send_jobs.write(' \n')
phases_line = 'phases='
for phase in phases_to_relax:
    phases_line = phases_line + '"' + phase + '" '
phases_line = phases_line + ' \n'
send_jobs.write(phases_line)
send_jobs.write(' \n')
send_jobs.write('for dir in "${phases[@]}"; do \n')
send_jobs.write('cd relax-phases/$dir \n')
send_jobs.write('sbatch run.sh \n')
send_jobs.write('cd ../.. \n')
send_jobs.write('done')
send_jobs.close()
os.system("chmod +x send_jobs.sh")