import os
import shutil

S_coef = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
Br_coef = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]

send_jobs = open('VCA_structures/send_jobs-bg.sh', "w")
send_jobs.write('#!/bin/bash \n')
send_jobs.write(' \n')

check_results = open('VCA_structures/check_results-bg.sh', "w")
check_results.write('#!/bin/bash \n')
check_results.write(' \n')

iteration = 1
for S_in in S_coef:
    for Br_in in Br_coef:
        # POSCAR
        source = 'VCA_structures/vca-' + str(iteration).zfill(3) + '/relaxation/CONTCAR'
        destination = 'VCA_structures/vca-' + str(iteration).zfill(3) + '/bandgap/POSCAR'
        shutil.copy(source, destination)

        # scripts
        send_jobs.write('cd ' + 'vca-' + str(iteration).zfill(3) + '/bandgap \n')
        send_jobs.write('sbatch run.sh \n')
        send_jobs.write('cd ../.. \n')
        send_jobs.write('\n')

        check_results.write('cd ' + 'vca-' + str(iteration).zfill(3) + '/bandgap \n')
        check_results.write('echo ' 'vca-' + str(iteration).zfill(3) + '\n')
        check_results.write('tail -n 1 OSZICAR \n')
        check_results.write('cd ../.. \n')
        check_results.write('\n')

        iteration = iteration + 1

os.system("chmod +x VCA_structures/send_jobs-bg.sh")
os.system("chmod +x VCA_structures/check_results-bg.sh")