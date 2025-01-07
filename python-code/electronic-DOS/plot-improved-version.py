import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

def get_pDOS(path_file):
    """
    Returns the s, p, and d density of electrons for the desired atom

    Inputs:
        path_file -> path to the vaspkit generated .dat file with the pDOS results 
    """

    energy = []
    s_dos = []
    p_dos = []
    d_dos = []

    with open(path_file, 'r') as f:
        f.readline()
        
        for line in f:
            energy.append(float(line.split()[0])+(0))

            s_dos.append(float(line.split()[1]))

            py = float(line.split()[5])
            pz = float(line.split()[5])
            px = float(line.split()[5])
            p_dos.append(py + pz + px)

            dxy = float(line.split()[5])
            dyz = float(line.split()[6])
            dz2 = float(line.split()[7])
            dxz = float(line.split()[8])
            dx2 = float(line.split()[9])
            d_dos.append(dxy + dyz + dz2 + dxz + dx2)

    return energy, s_dos, p_dos, d_dos


############## BaTiO3 ##############
path = '../../pDOS-candidates/BaTiO3/'
list_of_atoms = ['Ba', 'Ti', 'O']
list_of_colors = ['lightcoral', 'firebrick', 'red', 
                  'darkkhaki', 'forestgreen', 'lime',
                  'deepskyblue', 'royalblue', 'mediumpurple']

fig, axs = plt.subplots(2, 1, figsize=(4, 6))

axs[0].set_title('Non-distorted')
axs[1].set_title('Distorted')

axs[0].set_xlabel('$E-E_{F}$ (eV)')
axs[1].set_xlabel('$E-E_{F}$ (eV)')
axs[0].set_ylabel('eDOS (a.u.)')
axs[1].set_ylabel('eDOS (a.u.)')

axs[0].set_xlim(-4, 6)
axs[1].set_xlim(-4, 8)
axs[0].set_ylim(0, 3)
axs[1].set_ylim(0, 3)

it_color = 0
for atom in list_of_atoms:
    act_path = path + 'u0/PDOS_' + atom + '_SOC.dat'
    non_dist_energy, non_dist_s, non_dist_p, non_dist_d = get_pDOS(act_path)

    act_path = path + 'u4/PDOS_' + atom + '_SOC.dat'
    dist_energy, dist_s, dist_p, dist_d = get_pDOS(act_path)

    axs[0].plot(non_dist_energy, non_dist_s, color=list_of_colors[it_color], alpha=0.9, label=f'{atom}: $s$')
    axs[0].fill_between(non_dist_energy, non_dist_s, color=list_of_colors[it_color], alpha=0.3)
    axs[1].plot(dist_energy, dist_s, color=list_of_colors[it_color], alpha=0.9, label=f'{atom}: $s$')
    axs[1].fill_between(dist_energy, dist_s, color=list_of_colors[it_color], alpha=0.3)
    it_color = it_color + 1

    axs[0].plot(non_dist_energy, non_dist_p, color=list_of_colors[it_color], alpha=0.9, label=f'{atom}: $p$')
    axs[0].fill_between(non_dist_energy, non_dist_p, color=list_of_colors[it_color], alpha=0.3)
    axs[1].plot(dist_energy, dist_p, color=list_of_colors[it_color], alpha=0.9, label=f'{atom}: $p$')
    axs[1].fill_between(dist_energy, dist_p, color=list_of_colors[it_color], alpha=0.3)
    it_color = it_color + 1

    axs[0].plot(non_dist_energy, non_dist_d, color=list_of_colors[it_color], alpha=0.9, label=f'{atom}: $d$')
    axs[0].fill_between(non_dist_energy, non_dist_d, color=list_of_colors[it_color], alpha=0.3)
    axs[1].plot(dist_energy, dist_d, color=list_of_colors[it_color], alpha=0.9, label=f'{atom}: $d$')
    axs[1].fill_between(dist_energy, dist_d, color=list_of_colors[it_color], alpha=0.3)
    it_color = it_color + 1


major_locator = MultipleLocator(2) 
minor_locator = MultipleLocator(.4) 
axs[0].xaxis.set_major_locator(major_locator)
axs[0].xaxis.set_minor_locator(minor_locator)
major_locator = MultipleLocator(2) 
minor_locator = MultipleLocator(.4) 
axs[1].xaxis.set_major_locator(major_locator)
axs[1].xaxis.set_minor_locator(minor_locator)

major_locator = MultipleLocator(.5)
minor_locator = MultipleLocator(.1)
axs[0].yaxis.set_major_locator(major_locator)
axs[0].yaxis.set_minor_locator(minor_locator)
major_locator = MultipleLocator(.5)
minor_locator = MultipleLocator(.1)
axs[1].yaxis.set_major_locator(major_locator)
axs[1].yaxis.set_minor_locator(minor_locator)

axs[0].legend(frameon=False)

plt.tight_layout()
plt.savefig('PDOS_BaTiO3.pdf')
####################################


############## PbHfO3 ##############
path = '../../pDOS-candidates/PbHfO3/'
list_of_atoms = ['Pb', 'Hf', 'O']
list_of_colors = ['lightcoral', 'firebrick', 'red', 
                  'darkkhaki', 'forestgreen', 'lime',
                  'deepskyblue', 'royalblue', 'mediumpurple']

fig, axs = plt.subplots(2, 1, figsize=(4, 6))

axs[0].set_title('Non-distorted')
axs[1].set_title('Distorted')

axs[0].set_xlabel('$E-E_{F}$ (eV)')
axs[1].set_xlabel('$E-E_{F}$ (eV)')
axs[0].set_ylabel('eDOS (a.u.)')
axs[1].set_ylabel('eDOS (a.u.)')

axs[0].set_xlim(-4, 6)
axs[1].set_xlim(-4, 6)
axs[0].set_ylim(0, 3)
axs[1].set_ylim(0, 3)

it_color = 0
for atom in list_of_atoms:
    act_path = path + 'u0/PDOS_' + atom + '_SOC.dat'
    non_dist_energy, non_dist_s, non_dist_p, non_dist_d = get_pDOS(act_path)

    act_path = path + 'u4/PDOS_' + atom + '_SOC.dat'
    dist_energy, dist_s, dist_p, dist_d = get_pDOS(act_path)

    axs[0].plot(non_dist_energy, non_dist_s, color=list_of_colors[it_color], alpha=0.9, label=f'{atom}: $s$')
    axs[0].fill_between(non_dist_energy, non_dist_s, color=list_of_colors[it_color], alpha=0.3)
    axs[1].plot(dist_energy, dist_s, color=list_of_colors[it_color], alpha=0.9, label=f'{atom}: $s$')
    axs[1].fill_between(dist_energy, dist_s, color=list_of_colors[it_color], alpha=0.3)
    it_color = it_color + 1

    axs[0].plot(non_dist_energy, non_dist_p, color=list_of_colors[it_color], alpha=0.9, label=f'{atom}: $p$')
    axs[0].fill_between(non_dist_energy, non_dist_p, color=list_of_colors[it_color], alpha=0.3)
    axs[1].plot(dist_energy, dist_p, color=list_of_colors[it_color], alpha=0.9, label=f'{atom}: $p$')
    axs[1].fill_between(dist_energy, dist_p, color=list_of_colors[it_color], alpha=0.3)
    it_color = it_color + 1

    axs[0].plot(non_dist_energy, non_dist_d, color=list_of_colors[it_color], alpha=0.9, label=f'{atom}: $d$')
    axs[0].fill_between(non_dist_energy, non_dist_d, color=list_of_colors[it_color], alpha=0.3)
    axs[1].plot(dist_energy, dist_d, color=list_of_colors[it_color], alpha=0.9, label=f'{atom}: $d$')
    axs[1].fill_between(dist_energy, dist_d, color=list_of_colors[it_color], alpha=0.3)
    it_color = it_color + 1


major_locator = MultipleLocator(2) 
minor_locator = MultipleLocator(.4) 
axs[0].xaxis.set_major_locator(major_locator)
axs[0].xaxis.set_minor_locator(minor_locator)
major_locator = MultipleLocator(2) 
minor_locator = MultipleLocator(.4) 
axs[1].xaxis.set_major_locator(major_locator)
axs[1].xaxis.set_minor_locator(minor_locator)

major_locator = MultipleLocator(.5)
minor_locator = MultipleLocator(.1)
axs[0].yaxis.set_major_locator(major_locator)
axs[0].yaxis.set_minor_locator(minor_locator)
major_locator = MultipleLocator(.5)
minor_locator = MultipleLocator(.1)
axs[1].yaxis.set_major_locator(major_locator)
axs[1].yaxis.set_minor_locator(minor_locator)

axs[0].legend(frameon=False)

plt.tight_layout()
plt.savefig('PDOS_PbHfO3.pdf')
####################################