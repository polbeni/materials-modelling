# Pol Benítez Colominas, June 2025
# Universitat Politècnica de Catalunya

# Generate POSCAR files from snapshots of a MD simulation

import os
import shutil


################################ PARAMETERS ###############################
XDATCAR_path = 'XDATCAR'                                                # Path to the XDATCAR file

outputs_dir = 'outputs_file/'                                           # Path to dir where outputs are saved

total_configurations = 1200                                             # Total number of snapshots in the MD simulation
number_of_snapshots = 20                                                # Number of snapshots we want to consider
starting_point = 600                                                    # Initial snaphsot to consider
###########################################################################


def generate_POSCAR_from_snapshots(XDATCAR_path, starting_snap, snapshots_num, total_num, path_to_POSCARS):
    """
    Reads a XDATCAR file and save the snapshots as a POSCAR files

    Inputs:
        XDATCAR_path: path to the XDATCAR file
        starting_snap: first snapshot to consider
        snapshots_num: desired number of snapshots to save as POSCAR
        total_num: total number of snapshots in the XDATCAR file
        path_to_POSCARS: path of dir to save the POSCARS
    """

    same_lines = ['Snapshot from MLIP-MD\n']

    XDATCAR = open(XDATCAR_path, 'r')
    XDATCAR.readline()
    for _ in range(6):
        line = XDATCAR.readline()
        same_lines.append(line)
    XDATCAR.close()
    same_lines.append('Direct\n')

    num_atoms = 0
    for element in range(len(line.split())):
        num_atoms = num_atoms + int(line.split()[element])

    XDATCAR = open(XDATCAR_path, 'r')

    for _ in range(7):
        XDATCAR.readline()

    for snapshot in range(starting_snap):
        XDATCAR.readline()
        for _ in range(num_atoms):
            XDATCAR.readline()

    num_snapshot = 1
    for snapshot in range(total_num - starting_snap):
        XDATCAR.readline()

        if (snapshot % (int((total_num - starting_snap) / snapshots_num))) == 0:
            POSCAR = open(path_to_POSCARS + 'POSCAR-' + str(num_snapshot).zfill(3), 'w')
            for line in same_lines:
                POSCAR.write(line)

            for _ in range(num_atoms):
                line = XDATCAR.readline()
                POSCAR.write(line)

            POSCAR.close()
            num_snapshot = num_snapshot + 1
        else:
            for _ in range(num_atoms):
                XDATCAR.readline()

    XDATCAR.close()


if os.path.exists(outputs_dir):
    shutil.rmtree(outputs_dir)
os.mkdir(outputs_dir)

generate_POSCAR_from_snapshots(XDATCAR_path, starting_point, number_of_snapshots, total_configurations, outputs_dir)