# Pol Benítez Colominas, January 2024
# Universitat Politècnica de Catalunya

# Code to study the convergence after using parameters.py script
# For now the only parameters studied are the R2 of train and test set as well as the 
# loss function metric of some high symmetry points in the reciprocal-space, after the 
# harmonic phonon frequencies determination

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.gridspec import GridSpec
import yaml

from ase.io.vasp import read_vasp

from ase import Atoms
from ase.build import bulk
from hiphive import ForceConstantPotential

from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms

def loss_function(exact_freq, hiphive_freq):
    """
    Computes the loss function 

    Inputs:
        exact_freq: exact phonon frequencies at a given high symmetry point computed with Phonopy
        hiphive_freq: hiPhive determined phonon frequencies at a given high symmetry point
    """
    q_loss = 0

    for x in range(len(exact_freq)):
        q_loss = q_loss + (exact_freq[x] - hiphive_freq[x])**2

    return q_loss

def phonons_from_phonopy(file_path):
    """
    This function uses the band.yaml file generated in phonopy calculations to obtain the 
    phonon spectra data. This can be ploted with matplotlib

    Inputs:
        file_path: the path and name of the file

    Outputs:
        num_atoms: the number of atoms in the unit cell
        nqpoints: number of points in the reciprocal space
        npaths: the number of different paths in the reciprocal space
        segments_nqpoints: the number of points in the reciprocal space for each path
        points_labels: the labels of the high symmetry points
        phonons: a matrix array with the information about the phonons, it has the following structure
                    q-point   | branch 1 | branch 2 | branch 3 | ...
                    ------------------------------------------------
                    q-point 1 |   freq   |   freq   |   freq   | ...
                    q-point 2 |   freq   |   freq   |   freq   | ...
                    q-point 3 |   freq   |   freq   |   freq   | ...
                    ...       |   ...    |   ...    |   ...    | ...
    """
    with open(file_path, "r") as file:
        data = yaml.safe_load(file)

    num_atoms = data["natom"]
    nqpoints = data["nqpoint"]
    npaths = data["npath"]
    segments_nqpoints = data["segment_nqpoint"]
    points_labels = data["labels"]

    phonons = np.zeros((nqpoints, num_atoms*3 + 1))
    for x in range(nqpoints):
        phonons[x,0] = data["phonon"][x]["distance"]
        for y in range(num_atoms*3):
            phonons[x,y+1] = data["phonon"][x]["band"][y]["frequency"]

    return num_atoms, nqpoints, npaths, segments_nqpoints, points_labels, phonons

def get_band(q_start, q_stop, N):
    """ Return path between q_start and q_stop """
    return np.array([q_start + (q_stop-q_start)*i/(N-1) for i in range(N)])

### R2 convergence study

# create a file where storage the convergence results
results_r2 = open('convergence_r2.txt', 'w')
results_r2.write('cutoffs   #structures   R2_train   R2_test \n')

# define all the cutoffs of interest and number of structures
diff_cutoffs = [[4], [4.5], [5]]
diff_num_structures = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500]

# define an array with all the directories names
dir_results = ['1_050', '1_100', '1_150', '1_200', '1_250', '1_300', '1_350', '1_400', '1_450', '1_500',
               '2_050', '2_100', '2_150', '2_200', '2_250', '2_300', '2_350', '2_400', '2_450', '2_500',
               '3_050', '3_100', '3_150', '3_200', '3_250', '3_300', '3_350', '3_400', '3_450', '3_500']

num_iteration = 0

for cuts in diff_cutoffs:
    for struc in diff_num_structures:
        optm_file = open(dir_results[num_iteration] + '/log_files/4-optimizer', 'r')
        for x in range(11):
            optm_file.readline()

        actual_line = optm_file.readline()
        r2_train = actual_line.split()[2]
        actual_line = optm_file.readline()
        r2_test = actual_line.split()[2]

        results_r2.write(f'{cuts}   {struc}   {r2_train}   {r2_test} \n')

        optm_file.close()

        num_iteration = num_iteration + 1

results_r2.close()


### frequencies convergence study

dim = [2,2,2]  # dimension in phonopy calculation
Nq = 51  # number of q-points along each segment of the path through the BZ
mesh = [32, 32, 32]  # q-point mesh for MSD calculation

prim = read_vasp(file='POSCAR')
atoms_phonopy = PhonopyAtoms(symbols=prim.get_chemical_symbols(),
                             scaled_positions=prim.get_scaled_positions(),
                             cell=prim.cell)
phonopy = Phonopy(atoms_phonopy, supercell_matrix=dim*np.eye(3),
                  primitive_matrix=None)

G2X = get_band(np.array([0, 0, 0]), np.array([0.5, 0, 0]), Nq)
X2S = get_band(np.array([0.5, 0, 0]), np.array([0.5, 0.5, 0]), Nq)
S2Y = get_band(np.array([0.5, 0.5, 0]), np.array([0, 0.5, 0]), Nq)
Y2G = get_band(np.array([0, 0.5, 0]), np.array([0, 0, 0]), Nq)
bands = [G2X, X2S, S2Y, Y2G]

freqs_gamma_hiphive = [None]*45
freqs_x_hiphive = [None]*45
freqs_s_hiphive = [None]*45

num_iteration = 0
for cuts in diff_cutoffs:
    for struc in diff_num_structures:
        fcp = ForceConstantPotential.read(dir_results[num_iteration] + '/material.fcp')
        supercell = phonopy.get_supercell()
        supercell = Atoms(cell=supercell.cell, numbers=supercell.numbers, pbc=True,
                        scaled_positions=supercell.get_scaled_positions())
        fcs = fcp.get_force_constants(supercell)

        phonopy.set_force_constants(fcs.get_fc_array(order=2))
        phonopy.set_mesh(mesh, is_eigenvectors=True, is_mesh_symmetry=False)

        phonopy.set_band_structure(bands)
        qvecs, qnorms, freqs, _ = phonopy.get_band_structure()

        freqs_gamma_hiphive[num_iteration] = freqs[0][0,:]
        freqs_x_hiphive[num_iteration] = freqs[1][0,:]
        freqs_s_hiphive[num_iteration] = freqs[2][0,:]
        num_iteration = num_iteration + 1



# phonon frequencies from Phonopy
_, _, _, _, _, phonon_freq = phonons_from_phonopy('band.yaml')
freqs_gamma = [0]*36
freqs_x = [0]*36
freqs_s = [0]*36

for x in range(36):
    freqs_gamma[x] = phonon_freq[0,x]
    freqs_x[x] = phonon_freq[52,x]
    freqs_s[x] = phonon_freq[103,x]


results_loss = open('convergence_loss.txt', 'w')
results_loss.write('cutoffs   #structures   loss_Gamma   loss_X   loss_S \n')

num_iteration = 0
for cuts in diff_cutoffs:
    for struc in diff_num_structures:
        loss_value_gamma = loss_function(freqs_gamma, freqs_gamma_hiphive[num_iteration])
        loss_value_x = loss_function(freqs_x, freqs_x_hiphive[num_iteration])
        loss_value_s = loss_function(freqs_s, freqs_s_hiphive[num_iteration])

    

        results_loss.write(f'{cuts}   {struc}   {loss_value_gamma}   {loss_value_x}   {loss_value_s} \n')

        num_iteration = num_iteration + 1

results_loss.close()



### plot the results

fig, axs = plt.subplots(2, 1, figsize=(4,4))

axs[0].set_xlabel('# structures')
axs[0].set_ylabel('$R^2$ train')
axs[1].set_xlabel('# structures')
axs[1].set_ylabel('$R^2$ test')

num_structures = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500]

train_1 = []
train_2 = []
train_3 = []

test_1 = []
test_2 = []
test_3 = []

results_r2 = open('convergence_r2.txt', 'r')
results_r2.readline()

for x in range(len(num_structures)):
    actual_line = results_r2.readline()
    train_1.append(float(actual_line.split()[2]))
    test_1.append(float(actual_line.split()[3]))
    
for x in range(len(num_structures)):
    actual_line = results_r2.readline()
    train_2.append(float(actual_line.split()[2]))
    test_2.append(float(actual_line.split()[3]))

for x in range(len(num_structures)):
    actual_line = results_r2.readline()
    train_3.append(float(actual_line.split()[2]))
    test_3.append(float(actual_line.split()[3]))

results_r2.close()

axs[0].plot(num_structures, train_1, marker='o', linestyle='--', linewidth=0.5, 
            color='salmon', alpha=0.6, label='Cutoff: 4.0 Å')
axs[0].plot(num_structures, train_2, marker='o', linestyle='--', linewidth=0.5, 
            color='deepskyblue', alpha=0.6, label='Cutoff: 4.5 Å')
axs[0].plot(num_structures, train_3, marker='o', linestyle='--', linewidth=0.5, 
            color='lightgreen', alpha=0.6, label='Cutoff: 5.0 Å')

axs[1].plot(num_structures, test_1, marker='o', linestyle='--', linewidth=0.5, 
            color='salmon', alpha=0.6, label='Cutoff: 4.0 Å')
axs[1].plot(num_structures, test_2, marker='o', linestyle='--', linewidth=0.5, 
            color='deepskyblue', alpha=0.6, label='Cutoff: 4.5 Å')
axs[1].plot(num_structures, test_3, marker='o', linestyle='--', linewidth=0.5, 
            color='lightgreen', alpha=0.6, label='Cutoff: 5.0 Å')

axs[0].legend(fontsize=7)

major_locator = MultipleLocator(0.015) 
minor_locator = MultipleLocator(0.003) 
axs[0].yaxis.set_major_locator(major_locator)
axs[0].yaxis.set_minor_locator(minor_locator)
axs[1].yaxis.set_major_locator(major_locator)
axs[1].yaxis.set_minor_locator(minor_locator)

major_locator = MultipleLocator(100) 
minor_locator = MultipleLocator(50) 
axs[0].xaxis.set_major_locator(major_locator)
axs[0].xaxis.set_minor_locator(minor_locator)
axs[1].xaxis.set_major_locator(major_locator)
axs[1].xaxis.set_minor_locator(minor_locator)

axs[0].set_ylim(0.90, 0.96)
axs[1].set_ylim(0.92, 0.96)

axs[0].set_xlim(0, 550)
axs[1].set_xlim(0, 550)

plt.tight_layout()
plt.savefig('r2-convergence.pdf')



fig, axs = plt.subplots(3, 1, figsize=(6,6))

axs[0].set_xlabel('# structures')
axs[0].set_ylabel('Loss $\\Gamma$')
axs[1].set_xlabel('# structures')
axs[1].set_ylabel('Loss X')
axs[2].set_xlabel('# structures')
axs[2].set_ylabel('Loss S')

num_structures = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500]

gamma_1 = []
gamma_2 = []
gamma_3 = []

x_1 = []
x_2 = []
x_3 = []

s_1 = []
s_2 = []
s_3 = []

results_loss = open('convergence_loss.txt', 'r')
results_loss.readline()

for x in range(len(num_structures)):
    actual_line = results_loss.readline()
    gamma_1.append(float(actual_line.split()[2]))
    x_1.append(float(actual_line.split()[3]))
    s_1.append(float(actual_line.split()[4]))

for x in range(len(num_structures)):
    actual_line = results_loss.readline()
    gamma_2.append(float(actual_line.split()[2]))
    x_2.append(float(actual_line.split()[3]))
    s_2.append(float(actual_line.split()[4]))

for x in range(len(num_structures)):
    actual_line = results_loss.readline()
    gamma_3.append(float(actual_line.split()[2]))
    x_3.append(float(actual_line.split()[3]))
    s_3.append(float(actual_line.split()[4]))

results_loss.close()

axs[0].plot(num_structures, gamma_1, marker='o', linestyle='--', linewidth=0.5, 
            color='salmon', alpha=0.6, label='Cutoff: 4.0 Å')
axs[0].plot(num_structures, gamma_2, marker='o', linestyle='--', linewidth=0.5, 
            color='deepskyblue', alpha=0.6, label='Cutoff: 4.5 Å')
axs[0].plot(num_structures, gamma_3, marker='o', linestyle='--', linewidth=0.5, 
            color='lightgreen', alpha=0.6, label='Cutoff: 5.0 Å')

axs[1].plot(num_structures, x_1, marker='o', linestyle='--', linewidth=0.5, 
            color='salmon', alpha=0.6, label='Cutoff: 4.0 Å')
axs[1].plot(num_structures, x_2, marker='o', linestyle='--', linewidth=0.5, 
            color='deepskyblue', alpha=0.6, label='Cutoff: 4.5 Å')
axs[1].plot(num_structures, x_3, marker='o', linestyle='--', linewidth=0.5, 
            color='lightgreen', alpha=0.6, label='Cutoff: 5.0 Å')

axs[2].plot(num_structures, s_1, marker='o', linestyle='--', linewidth=0.5, 
            color='salmon', alpha=0.6, label='Cutoff: 4.0 Å')
axs[2].plot(num_structures, s_2, marker='o', linestyle='--', linewidth=0.5, 
            color='deepskyblue', alpha=0.6, label='Cutoff: 4.5 Å')
axs[2].plot(num_structures, s_3, marker='o', linestyle='--', linewidth=0.5, 
            color='lightgreen', alpha=0.6, label='Cutoff: 5.0 Å')

axs[0].legend(fontsize=7)

major_locator = MultipleLocator(100) 
minor_locator = MultipleLocator(50) 
axs[0].xaxis.set_major_locator(major_locator)
axs[0].xaxis.set_minor_locator(minor_locator)
axs[1].xaxis.set_major_locator(major_locator)
axs[1].xaxis.set_minor_locator(minor_locator)
axs[2].xaxis.set_major_locator(major_locator)
axs[2].xaxis.set_minor_locator(minor_locator)

plt.tight_layout()
plt.savefig('loss-convergence.pdf')