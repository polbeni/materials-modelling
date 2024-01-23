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


#hiphive

dim = [3,3,2]  # dimension in phonopy calculation
Nq = 51  # number of q-points along each segment of the path through the BZ
mesh = [32, 32, 32]  # q-point mesh for MSD calculation

prim = read_vasp(file='POSCAR')
atoms_phonopy = PhonopyAtoms(symbols=prim.get_chemical_symbols(),
                             scaled_positions=prim.get_scaled_positions(),
                             cell=prim.cell)
phonopy = Phonopy(atoms_phonopy, supercell_matrix=dim*np.eye(3),
                  primitive_matrix=None)

fcp = ForceConstantPotential.read('material.fcp')
supercell = phonopy.get_supercell()
supercell = Atoms(cell=supercell.cell, numbers=supercell.numbers, pbc=True,
                  scaled_positions=supercell.get_scaled_positions())
fcs = fcp.get_force_constants(supercell)

phonopy.set_force_constants(fcs.get_fc_array(order=2))
phonopy.set_mesh(mesh, is_eigenvectors=True, is_mesh_symmetry=False)


def get_band(q_start, q_stop, N):
    """ Return path between q_start and q_stop """
    return np.array([q_start + (q_stop-q_start)*i/(N-1) for i in range(N)])


G2X = get_band(np.array([0, 0, 0]), np.array([0.5, 0, 0]), Nq)
X2S = get_band(np.array([0.5, 0, 0]), np.array([0.5, 0.5, 0]), Nq)
S2Y = get_band(np.array([0.5, 0.5, 0]), np.array([0, 0.5, 0]), Nq)
Y2G = get_band(np.array([0, 0.5, 0]), np.array([0, 0, 0]), Nq)
bands = [G2X, X2S, S2Y, Y2G]

# get phonon dispersion
phonopy.set_band_structure(bands)
qvecs, qnorms, freqs, _ = phonopy.get_band_structure()


fig = plt.figure()
plt.title('ZrO$_2$, $P4_2/nmc$')

kpts = [0.0, qnorms[0][-1], qnorms[1][-1], qnorms[2][-1], qnorms[3][-1]]
kpts_labels = ['$\\Gamma$', 'X', 'S', 'Y', '$\\Gamma$']

plt.axvline(x=kpts[1], color='k', linewidth=0.9)
plt.axvline(x=kpts[2], color='k', linewidth=0.9)
plt.axvline(x=kpts[3], color='k', linewidth=0.9)

plt.plot(qnorms[0], freqs[0][:,0], color='darkorange', alpha=0.9, linewidth=.8, label='Phonopy+hiPhive')
for q, freq, in zip(qnorms, freqs):
    plt.plot(q, freq, color='darkorange', alpha=0.9, linewidth=.8)

plt.ylabel('Freq (THz)', fontsize=14.0)
plt.xticks(kpts, kpts_labels, fontsize=14.0)
plt.xlim([0.0, qnorms[-1][-1]])
#plt.ylim([0.0, 12.0])


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

phase = 'phase-1'
file_path = 'band.yaml'

num_atoms, nqpoints, npaths, segments_nqpoints, points_labels, phonons = phonons_from_phonopy(file_path)

vertical_lines = ['False']*(npaths-1)
for x in range(npaths-1):
    if points_labels[x][1] != points_labels[x+1][0]:
        vertical_lines[x] = True
        
x_labels = ['point']*(npaths+1)
for x in range(npaths):
    if x == 0:
        x_labels[0] = points_labels[0][0]
    if  x+1 < npaths:
        if vertical_lines[x] == True:
            x_labels[x+1] = points_labels[x][1] + '|' + points_labels[x+1][0]
            if phase == 'phase-2':
                x_labels[x+1] = points_labels[x][1] + '\n' + points_labels[x+1][0]
        else: 
            x_labels[x+1] = points_labels[x][1]
    else:
        x_labels[x+1] = points_labels[x][1]


plt.axhline(0, color='black', linestyle='--', linewidth=1)

plt.plot(phonons[:,0], phonons[:,1], color='black', alpha=1, linewidth=.9, label='Phonopy')
for x in range(num_atoms*3):
    plt.plot(phonons[:,0], phonons[:,x+1], color='black', alpha=1, linewidth=.9)

plt.ylim(-3,27)


plt.legend(loc='upper right')
plt.tight_layout()
plt.savefig('phonon_dispersion.pdf')