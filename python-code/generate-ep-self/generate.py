# Pol Benítez Colominas, November 2024
# University of Cambridge and Universitat Politècnica de Catalunya

# Generates the files to calculate the electron-phonon matrix elements for self-scattering

import os
import shutil

import numpy as np

num_branch = [1, 2, 3, 4, 5]
num_point = [1, 2, 3, 4, 5, 6, 7]
dist_point = [0.0075, 0.015, 0.0225, 0.03, 0.0375, 0.045, 0.0525]

######################## SIGMA ########################
path = 'Sigma/'

send_jobs = open(path + 'send_jobs.sh', "w")
send_jobs.write('#!/bin/bash \n')
send_jobs.write(' \n')

check_results = open(path + 'check_results.sh', "w")
check_results.write('#!/bin/bash \n')
check_results.write(' \n')

for branch in num_branch:
    for point in num_point:
        os.mkdir(path + str(branch).zfill(1) + '_' + str(point).zfill(1))

        # copy scf.in file
        source = 'save/scf.in' 
        destination = path + str(branch).zfill(1) + '_' + str(point).zfill(1) + '/scf.in'
        shutil.copy(source, destination)

        # copy run.sh file
        source = 'save/run.sh' 
        destination = path + str(branch).zfill(1) + '_' + str(point).zfill(1) + '/run.sh'
        shutil.copy(source, destination)

        # copy pseudo files
        os.mkdir(path + str(branch).zfill(1) + '_' + str(point).zfill(1) + '/pseudo')

        source = 'save/pseudo/Ti-nc-sr-04_pbesol_stringent.upf' 
        destination = path + str(branch).zfill(1) + '_' + str(point).zfill(1) + '/pseudo/Ti-nc-sr-04_pbesol_stringent.upf'
        shutil.copy(source, destination)

        source = 'save/pseudo/S-nc-sr-04_pbesol_stringent.upf' 
        destination = path + str(branch).zfill(1) + '_' + str(point).zfill(1) + '/pseudo/S-nc-sr-04_pbesol_stringent.upf'
        shutil.copy(source, destination)

        # generate ph.in file
        ph_file = open(path + str(branch).zfill(1) + '_' + str(point).zfill(1) + '/ph.in', 'w')
        ph_file.write('\n')
        ph_file.write('&inputph\n')
        ph_file.write('\n')
        ph_file.write('tr2_ph=1.0d-15\n')
        ph_file.write('amass(1) = 47.867\n')
        ph_file.write('amass(2) = 32.065\n')
        ph_file.write("prefix='TiS2'\n")
        ph_file.write("outdir='tmp_dir'\n")
        ph_file.write("fildyn='TiS2_elecphonon.dyn'\n")
        ph_file.write("fildvscf = 'dvscf'\n")
        ph_file.write("electron_phonon='prt'\n")
        ph_file.write('kx = 0.0000000\n')
        ph_file.write('ky = 0.4811252\n')
        ph_file.write('kz = 0.0000000\n')
        ph_file.write('/\n')
        if branch == 1:
            ph_file.write(f'0.0 {dist_point[point - 1]:.6f} 0.0')
        elif branch == 2:
            ph_file.write(f'{(dist_point[point - 1] / np.sqrt(2)):.6f} {(dist_point[point - 1] / np.sqrt(2)):.6f} 0.0')
        elif branch == 3:
            ph_file.write(f'{dist_point[point - 1]:.6f} 0.0 0.0')
        elif branch == 4:
            ph_file.write(f'{(dist_point[point - 1] / np.sqrt(2)):.6f} -{(dist_point[point - 1] / np.sqrt(2)):.6f} 0.0')
        elif branch == 5:
            ph_file.write(f'0.0 -{dist_point[point - 1]:.6f} 0.0')
        ph_file.close()

        send_jobs.write('cd ' + str(branch).zfill(1) + '_' + str(point).zfill(1) + '\n')
        send_jobs.write('sbatch run.sh \n')
        send_jobs.write('cd .. \n')
        send_jobs.write('\n')

        check_results.write('cd ' + str(branch).zfill(1) + '_' + str(point).zfill(1) + '\n')
        check_results.write('tail ph.out \n')
        check_results.write('cd .. \n')
        check_results.write('\n')

send_jobs.close()
os.system("chmod +x " + path + "send_jobs.sh")

check_results.close()
os.system("chmod +x " + path + "check_results.sh")
#######################################################


######################## LAMBDA #######################
path = 'Lambda/'

def rotation(coord_x, coord_y, alpha):
    """
    It rotates a point in cartesian coordinates (in units 2pi/alat)

    Inputs:
        coord_x -> coordinate x
        coord_y -> coordinate y
        alpha -> degree of rotation (in rad)
    """

    new_x = (np.cos(alpha) * coord_x) + (-np.sin(alpha) * coord_y)
    new_y = (np.sin(alpha) * coord_x) + (np.cos(alpha) * coord_y)

    return new_x, new_y

send_jobs = open(path + 'send_jobs.sh', "w")
send_jobs.write('#!/bin/bash \n')
send_jobs.write(' \n')

check_results = open(path + 'check_results.sh', "w")
check_results.write('#!/bin/bash \n')
check_results.write(' \n')

for branch in num_branch:
    for point in num_point:
        os.mkdir(path + str(branch).zfill(1) + '_' + str(point).zfill(1))

        # copy scf.in file
        source = 'save/scf.in' 
        destination = path + str(branch).zfill(1) + '_' + str(point).zfill(1) + '/scf.in'
        shutil.copy(source, destination)

        # copy run.sh file
        source = 'save/run.sh' 
        destination = path + str(branch).zfill(1) + '_' + str(point).zfill(1) + '/run.sh'
        shutil.copy(source, destination)

        # copy pseudo files
        os.mkdir(path + str(branch).zfill(1) + '_' + str(point).zfill(1) + '/pseudo')

        source = 'save/pseudo/Ti-nc-sr-04_pbesol_stringent.upf' 
        destination = path + str(branch).zfill(1) + '_' + str(point).zfill(1) + '/pseudo/Ti-nc-sr-04_pbesol_stringent.upf'
        shutil.copy(source, destination)

        source = 'save/pseudo/S-nc-sr-04_pbesol_stringent.upf' 
        destination = path + str(branch).zfill(1) + '_' + str(point).zfill(1) + '/pseudo/S-nc-sr-04_pbesol_stringent.upf'
        shutil.copy(source, destination)

        # generate ph.in file
        ph_file = open(path + str(branch).zfill(1) + '_' + str(point).zfill(1) + '/ph.in', 'w')
        ph_file.write('\n')
        ph_file.write('&inputph\n')
        ph_file.write('\n')
        ph_file.write('tr2_ph=1.0d-15\n')
        ph_file.write('amass(1) = 47.867\n')
        ph_file.write('amass(2) = 32.065\n')
        ph_file.write("prefix='TiS2'\n")
        ph_file.write("outdir='tmp_dir'\n")
        ph_file.write("fildyn='TiS2_elecphonon.dyn'\n")
        ph_file.write("fildvscf = 'dvscf'\n")
        ph_file.write("electron_phonon='prt'\n")
        ph_file.write('kx = 0.0833333\n')
        ph_file.write('ky = 0.1443376\n')
        ph_file.write('kz = 0.0000000\n')
        ph_file.write('/\n')
        if branch == 1:
            val_x, val_y = rotation(0, dist_point[point - 1], -0.523599)
            ph_file.write(f'{val_x:.6f} {val_y:.6f} 0.0')
        elif branch == 2:
            val_x, val_y = rotation((dist_point[point - 1] / np.sqrt(2)), (dist_point[point - 1] / np.sqrt(2)), -0.523599)
            ph_file.write(f'{val_x:.6f} {val_y:.6f} 0.0')
        elif branch == 3:
            val_x, val_y = rotation(dist_point[point - 1], 0, -0.523599)
            ph_file.write(f'{val_x:.6f} {val_y:.6f} 0.0')
        elif branch == 4:
            val_x, val_y = rotation((dist_point[point - 1] / np.sqrt(2)), -(dist_point[point - 1] / np.sqrt(2)), -0.523599)
            ph_file.write(f'{val_x:.6f} {val_y:.6f} 0.0')
        elif branch == 5:
            val_x, val_y = rotation(0, -dist_point[point - 1], -0.523599)
            ph_file.write(f'{val_x:.6f} {val_y:.6f} 0.0')
        ph_file.close()

        send_jobs.write('cd ' + str(branch).zfill(1) + '_' + str(point).zfill(1) + '\n')
        send_jobs.write('sbatch run.sh \n')
        send_jobs.write('cd .. \n')
        send_jobs.write('\n')

        check_results.write('cd ' + str(branch).zfill(1) + '_' + str(point).zfill(1) + '\n')
        check_results.write('tail ph.out \n')
        check_results.write('cd .. \n')
        check_results.write('\n')

send_jobs.close()
os.system("chmod +x " + path + "send_jobs.sh")

check_results.close()
os.system("chmod +x " + path + "check_results.sh")
#######################################################