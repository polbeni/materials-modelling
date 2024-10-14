import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.gridspec import GridSpec
from scipy.interpolate import interp1d

plt.figure()
fig, axs = plt.subplots(4, 2, figsize=(4,8))

gs = GridSpec(4, 2, width_ratios=[1, 0.3])


axs[0,0] = plt.subplot(gs[0,0])
axs[0,1] = plt.subplot(gs[0,1])
axs[1,0] = plt.subplot(gs[1,0])
axs[1,1] = plt.subplot(gs[1,1])
axs[2,0] = plt.subplot(gs[2,0])
axs[2,1] = plt.subplot(gs[2,1])
axs[3,0] = plt.subplot(gs[3,0])
axs[3,1] = plt.subplot(gs[3,1])

major_locator = MultipleLocator(5) 
minor_locator = MultipleLocator(1) 
axs[0,0].yaxis.set_major_locator(major_locator)
axs[0,0].yaxis.set_minor_locator(minor_locator)
axs[1,0].yaxis.set_major_locator(major_locator)
axs[1,0].yaxis.set_minor_locator(minor_locator)
axs[2,0].yaxis.set_major_locator(major_locator)
axs[2,0].yaxis.set_minor_locator(minor_locator)
axs[3,0].yaxis.set_major_locator(major_locator)
axs[3,0].yaxis.set_minor_locator(minor_locator)

axs[0,1].yaxis.set_major_locator(major_locator)
axs[0,1].yaxis.set_minor_locator(minor_locator)
axs[1,1].yaxis.set_major_locator(major_locator)
axs[1,1].yaxis.set_minor_locator(minor_locator)
axs[2,1].yaxis.set_major_locator(major_locator)
axs[2,1].yaxis.set_minor_locator(minor_locator)
axs[3,1].yaxis.set_major_locator(major_locator)
axs[3,1].yaxis.set_minor_locator(minor_locator)



##### Ag3SBr, Pm3m

material = 'Ag3SBr-Pm3m'

axs[0,0].set_title('Ag$_3$SBr, $Pm\overline{3}m$', fontsize=12)
axs[0,0].set_ylabel('Energy (eV)')
axs[0,0].set_xlim(0, 1)
axs[0,0].set_ylim(-7.5, 12.5)


kpoint1 = np.zeros(192)
kpoint2 = np.zeros(192)
kpoint3 = np.zeros(192)
kpoint4 = np.zeros(192)
kpoint5 = np.zeros(192)
kpoint6 = np.zeros(192)
kpoint7 = np.zeros(192)
kpoint8 = np.zeros(192)
kpoint9 = np.zeros(192)
kpoint10 = np.zeros(192)
kpoint11 = np.zeros(192)
kpoint12 = np.zeros(192)
kpoint13 = np.zeros(192)
kpoint14 = np.zeros(192)
kpoint15 = np.zeros(192)
kpoint16 = np.zeros(192)
kpoint17 = np.zeros(192)
kpoint18 = np.zeros(192)
kpoint19 = np.zeros(192)
kpoint20 = np.zeros(192)
kpoints = [kpoint1, kpoint2, kpoint3, kpoint4, kpoint5, kpoint6, kpoint7, kpoint8, kpoint9, kpoint10,
           kpoint11, kpoint12, kpoint13, kpoint14, kpoint15, kpoint16, kpoint17, kpoint18, kpoint19, kpoint20]

path = 'data/' + material + '/OUTCAR'

OUTCAR_file = open(path, "r")
for x in range(1768):
    OUTCAR_file.readline()

for kpoint in kpoints:
    for x in range(192):
        actual_line = OUTCAR_file.readline()
        energy_value = actual_line.split()[1]

        kpoint[x] = energy_value

    for x in range(3):
        OUTCAR_file.readline()

OUTCAR_file.close()

band_matrix = np.zeros((192, 16))
y = 0

points_ordered = [1, 2, 3, 4, 7, 9, 10, 8, 5, 1, 11, 17, 20, 18, 13, 4]

for kpoint in points_ordered:
    for x in range(192):
        band_matrix[x,y] = kpoints[kpoint-1][x]
    y = y + 1




kspace1 = np.linspace(0, 1, 16)
kspace2 = np.linspace(1, 1.2, 4)

interv1 = np.linspace(0, 1, 160)
interv2 = np.linspace(1, 1.2, 16)
k_space_interp = np.linspace(0, 1, 160)

band_matrix_interp = np.zeros((192, 160))

for x in range(192):
    band_matrix_interp[x,:] = interp1d(kspace1, band_matrix[x,:],
                                  kind='quadratic', fill_value='extrapolate')(interv1)


x_labels = ['$\\Gamma$', 'X', 'M', '$\\Gamma$', 'R', 'X']
x_ticks = [0, 0.2, 0.4, 0.6, 0.8, 1]
for x in range(5):
    axs[0,0].axvline(x_ticks[x], color='black', linestyle='--', linewidth=.8)
axs[0,0].set_xticks(ticks=x_ticks, labels=x_labels)

for x in range(192):
    if x <= 22:
        color_type = 'coral'
    else:
        color_type = 'slateblue'
    axs[0,0].plot(interv1, band_matrix_interp[x,:]-np.max(band_matrix[0:22,:]), color=color_type)

for x in range(22):
    for y in range(16):
        if band_matrix[x,y] == np.max(band_matrix[0:22,:]):
            ind_max = y
axs[0,0].plot(kspace1[ind_max], np.max(band_matrix[0:22,:])-np.max(band_matrix[0:22,:]), 'o', markersize=8, color='darkred')

for x in range(23,192):
    for y in range(16):
        if band_matrix[x,y] == np.min(band_matrix[23:-1,:]):
            ind_min = y
axs[0,0].plot(kspace1[ind_min], np.min(band_matrix[23:-1,:])-np.max(band_matrix[0:22,:]), 'o', markersize=8, color='navy')

path = 'data/' + material + '/TDOS.dat'
total_dos = np.loadtxt(path)
path = 'data/' + material + '/PDOS_Ag.dat'
Ag_dos = np.loadtxt(path)
path = 'data/' + material + '/PDOS_S.dat'
S_dos = np.loadtxt(path)
path = 'data/' + material + '/PDOS_Br.dat'
Br_dos = np.loadtxt(path)

axs[0,1].plot(total_dos[:,-1], total_dos[:,0], color='black', linewidth=.9)
axs[0,1].fill_between(total_dos[:,-1], total_dos[:,0], color='grey', edgecolor='None', alpha=0.5)
axs[0,1].plot(Ag_dos[:,-1], Ag_dos[:,0], color='cornflowerblue', linewidth=.9)
axs[0,1].fill_between(Ag_dos[:,-1], Ag_dos[:,0], color='cornflowerblue', edgecolor='None', alpha=0.5)
axs[0,1].plot(S_dos[:,-1], S_dos[:,0], color='lime', linewidth=.9)
axs[0,1].fill_between(S_dos[:,-1], S_dos[:,0], color='lime', edgecolor='None', alpha=0.5)
axs[0,1].plot(Br_dos[:,-1], Br_dos[:,0], color='violet', linewidth=.9)
axs[0,1].fill_between(Br_dos[:,-1], Br_dos[:,0], color='violet', edgecolor='None', alpha=0.5)
axs[0,1].set_xlim(0,6)
axs[0,1].set_ylim(-7.5, 12.5)
axs[0,1].set_yticklabels([])

major_locator = MultipleLocator(3) 
minor_locator = MultipleLocator(.6)  

axs[0,1].xaxis.set_major_locator(major_locator)
axs[0,1].xaxis.set_minor_locator(minor_locator)



##### Ag3SBr, gs

material = 'Ag3SBr-gs'

axs[1,0].set_title('Ag$_3$SBr, $P2_13$', fontsize=12)
axs[1,0].set_ylabel('Energy (eV)')
axs[1,0].set_xlim(0, 1)
axs[1,0].set_ylim(-5, 7.5)


kpoint1 = np.zeros(192)
kpoint2 = np.zeros(192)
kpoint3 = np.zeros(192)
kpoint4 = np.zeros(192)
kpoint5 = np.zeros(192)
kpoint6 = np.zeros(192)
kpoint7 = np.zeros(192)
kpoint8 = np.zeros(192)
kpoint9 = np.zeros(192)
kpoint10 = np.zeros(192)
kpoint11 = np.zeros(192)
kpoints = [kpoint1, kpoint2, kpoint3, kpoint4, kpoint5, kpoint6, kpoint7, kpoint8, kpoint9, kpoint10,
           kpoint11]

path = 'data/' + material + '/OUTCAR'

OUTCAR_file = open(path, "r")
for x in range(1699):
    OUTCAR_file.readline()

for kpoint in kpoints:
    for x in range(192):
        actual_line = OUTCAR_file.readline()
        energy_value = actual_line.split()[1]

        kpoint[x] = energy_value

    for x in range(3):
        OUTCAR_file.readline()

OUTCAR_file.close()

band_matrix = np.zeros((192, 11))
y = 0

points_ordered = [1, 2, 3, 5, 7, 4, 1, 8, 11, 9, 3]

for kpoint in points_ordered:
    for x in range(192):
        band_matrix[x,y] = kpoints[kpoint-1][x]
    y = y + 1




kspace1 = np.linspace(0, 1, 11)
interv1 = np.linspace(0, 1, 110)
k_space_interp = np.linspace(0, 1, 110)

band_matrix_interp = np.zeros((192, 110))

for x in range(192):
    band_matrix_interp[x,:] = interp1d(kspace1, band_matrix[x,:],
                                  kind='quadratic', fill_value='extrapolate')(interv1)


x_labels = ['$\\Gamma$', 'X', 'M', '$\\Gamma$', 'R', 'X']
x_ticks = [0, 0.2, 0.4, 0.6, 0.8, 1]
for x in range(5):
    axs[1,0].axvline(x_ticks[x], color='black', linestyle='--', linewidth=.8)
axs[1,0].set_xticks(ticks=x_ticks, labels=x_labels)

for x in range(192):
    if x <= 91:
        color_type = 'coral'
    else:
        color_type = 'slateblue'
    axs[1,0].plot(interv1, band_matrix_interp[x,:]-np.max(band_matrix[0:92,:]), color=color_type)

for x in range(92):
    for y in range(11):
        if band_matrix[x,y] == np.max(band_matrix[0:92,:]):
            ind_max = y
axs[1,0].plot(kspace1[ind_max], np.max(band_matrix[0:92,:])-np.max(band_matrix[0:92,:]), 'o', markersize=8, color='darkred')

for x in range(92,192):
    for y in range(11):
        if band_matrix[x,y] == np.min(band_matrix[92:-1,:]):
            ind_min = y
axs[1,0].plot(kspace1[ind_min], np.min(band_matrix[92:-1,:])-np.max(band_matrix[0:92,:]), 'o', markersize=8, color='navy')


path = 'data/' + material + '/TDOS.dat'
total_dos = np.loadtxt(path)
path = 'data/' + material + '/PDOS_Ag.dat'
Ag_dos = np.loadtxt(path)
path = 'data/' + material + '/PDOS_S.dat'
S_dos = np.loadtxt(path)
path = 'data/' + material + '/PDOS_Br.dat'
Br_dos = np.loadtxt(path)

axs[1,1].plot(total_dos[:,-1], total_dos[:,0], color='black', linewidth=.9)
axs[1,1].fill_between(total_dos[:,-1], total_dos[:,0], color='grey', edgecolor='None', alpha=0.5)
axs[1,1].plot(Ag_dos[:,-1], Ag_dos[:,0], color='cornflowerblue', linewidth=.9)
axs[1,1].fill_between(Ag_dos[:,-1], Ag_dos[:,0], color='cornflowerblue', edgecolor='None', alpha=0.5)
axs[1,1].plot(S_dos[:,-1], S_dos[:,0], color='lime', linewidth=.9)
axs[1,1].fill_between(S_dos[:,-1], S_dos[:,0], color='lime', edgecolor='None', alpha=0.5)
axs[1,1].plot(Br_dos[:,-1], Br_dos[:,0], color='violet', linewidth=.9)
axs[1,1].fill_between(Br_dos[:,-1], Br_dos[:,0], color='violet', edgecolor='None', alpha=0.5)
axs[1,1].set_xlim(0,45)
axs[1,1].set_ylim(-5,7.5)
axs[1,1].set_yticklabels([])

major_locator = MultipleLocator(20) 
minor_locator = MultipleLocator(4)  

axs[1,1].xaxis.set_major_locator(major_locator)
axs[1,1].xaxis.set_minor_locator(minor_locator)








##### Ag3SI, Pm3m

material = 'Ag3SI-Pm3m'

axs[2,0].set_title('Ag$_3$SI, $Pm\overline{3}m$', fontsize=12)
axs[2,0].set_ylabel('Energy (eV)')
axs[2,0].set_xlim(0, 1)
axs[2,0].set_ylim(-5, 7.5)


kpoint1 = np.zeros(192)
kpoint2 = np.zeros(192)
kpoint3 = np.zeros(192)
kpoint4 = np.zeros(192)
kpoint5 = np.zeros(192)
kpoint6 = np.zeros(192)
kpoint7 = np.zeros(192)
kpoint8 = np.zeros(192)
kpoint9 = np.zeros(192)
kpoint10 = np.zeros(192)
kpoint11 = np.zeros(192)
kpoint12 = np.zeros(192)
kpoint13 = np.zeros(192)
kpoint14 = np.zeros(192)
kpoint15 = np.zeros(192)
kpoint16 = np.zeros(192)
kpoint17 = np.zeros(192)
kpoint18 = np.zeros(192)
kpoint19 = np.zeros(192)
kpoint20 = np.zeros(192)
kpoints = [kpoint1, kpoint2, kpoint3, kpoint4, kpoint5, kpoint6, kpoint7, kpoint8, kpoint9, kpoint10,
           kpoint11, kpoint12, kpoint13, kpoint14, kpoint15, kpoint16, kpoint17, kpoint18, kpoint19, kpoint20]

path = 'data/' + material + '/OUTCAR'

OUTCAR_file = open(path, "r")
for x in range(1772):
    OUTCAR_file.readline()

for kpoint in kpoints:
    for x in range(192):
        actual_line = OUTCAR_file.readline()
        energy_value = actual_line.split()[1]

        kpoint[x] = energy_value

    for x in range(3):
        OUTCAR_file.readline()

OUTCAR_file.close()

band_matrix = np.zeros((192, 16))
y = 0

points_ordered = [1, 2, 3, 4, 7, 9, 10, 8, 5, 1, 11, 17, 20, 18, 13, 4]

for kpoint in points_ordered:
    for x in range(192):
        band_matrix[x,y] = kpoints[kpoint-1][x]
    y = y + 1




kspace1 = np.linspace(0, 1, 16)
kspace2 = np.linspace(1, 1.2, 4)

interv1 = np.linspace(0, 1, 160)
interv2 = np.linspace(1, 1.2, 16)
k_space_interp = np.linspace(0, 1, 160)

band_matrix_interp = np.zeros((192, 160))

for x in range(192):
    band_matrix_interp[x,:] = interp1d(kspace1, band_matrix[x,:],
                                  kind='quadratic', fill_value='extrapolate')(interv1)


x_labels = ['$\\Gamma$', 'X', 'M', '$\\Gamma$', 'R', 'X']
x_ticks = [0, 0.2, 0.4, 0.6, 0.8, 1]
for x in range(5):
    axs[2,0].axvline(x_ticks[x], color='black', linestyle='--', linewidth=.8)
axs[2,0].set_xticks(ticks=x_ticks, labels=x_labels)

for x in range(192):
    if x <= 22:
        color_type = 'coral'
    else:
        color_type = 'slateblue'
    axs[2,0].plot(interv1, band_matrix_interp[x,:]-np.max(band_matrix[0:22,:]), color=color_type)

for x in range(22):
    for y in range(16):
        if band_matrix[x,y] == np.max(band_matrix[0:22,:]):
            ind_max = y
axs[2,0].plot(kspace1[ind_max], np.max(band_matrix[0:22,:])-np.max(band_matrix[0:22,:]), 'o', markersize=8, color='darkred')

for x in range(23,192):
    for y in range(16):
        if band_matrix[x,y] == np.min(band_matrix[23:-1,:]):
            ind_min = y
axs[2,0].plot(kspace1[ind_min], np.min(band_matrix[23:-1,:])-np.max(band_matrix[0:22,:]), 'o', markersize=8, color='navy')

path = 'data/' + material + '/TDOS.dat'
total_dos = np.loadtxt(path)
path = 'data/' + material + '/PDOS_Ag.dat'
Ag_dos = np.loadtxt(path)
path = 'data/' + material + '/PDOS_S.dat'
S_dos = np.loadtxt(path)
path = 'data/' + material + '/PDOS_I.dat'
Br_dos = np.loadtxt(path)

axs[2,1].plot(total_dos[:,-1], total_dos[:,0], color='black', linewidth=.9)
axs[2,1].fill_between(total_dos[:,-1], total_dos[:,0], color='grey', edgecolor='None', alpha=0.5)
axs[2,1].plot(Ag_dos[:,-1], Ag_dos[:,0], color='cornflowerblue', linewidth=.9)
axs[2,1].fill_between(Ag_dos[:,-1], Ag_dos[:,0], color='cornflowerblue', edgecolor='None', alpha=0.5)
axs[2,1].plot(S_dos[:,-1], S_dos[:,0], color='lime', linewidth=.9)
axs[2,1].fill_between(S_dos[:,-1], S_dos[:,0], color='lime', edgecolor='None', alpha=0.5)
axs[2,1].plot(Br_dos[:,-1], Br_dos[:,0], color='violet', linewidth=.9)
axs[2,1].fill_between(Br_dos[:,-1], Br_dos[:,0], color='violet', edgecolor='None', alpha=0.5)
axs[2,1].set_xlim(0,9)
axs[2,1].set_ylim(-5,7.5)
axs[2,1].set_yticklabels([])

major_locator = MultipleLocator(4) 
minor_locator = MultipleLocator(.8)  

axs[2,1].xaxis.set_major_locator(major_locator)
axs[2,1].xaxis.set_minor_locator(minor_locator)



##### Ag3SI, gs

material = 'Ag3SI-gs'

axs[3,0].set_title('Ag$_3$SI, $P2_13$', fontsize=12)
axs[3,0].set_ylabel('Energy (eV)')
axs[3,0].set_xlim(0, 1)
axs[3,0].set_ylim(-5, 7.5)


kpoint1 = np.zeros(192)
kpoint2 = np.zeros(192)
kpoint3 = np.zeros(192)
kpoint4 = np.zeros(192)
kpoint5 = np.zeros(192)
kpoint6 = np.zeros(192)
kpoint7 = np.zeros(192)
kpoint8 = np.zeros(192)
kpoint9 = np.zeros(192)
kpoint10 = np.zeros(192)
kpoint11 = np.zeros(192)
kpoints = [kpoint1, kpoint2, kpoint3, kpoint4, kpoint5, kpoint6, kpoint7, kpoint8, kpoint9, kpoint10,
           kpoint11]

path = 'data/' + material + '/OUTCAR'

OUTCAR_file = open(path, "r")
for x in range(1714):
    OUTCAR_file.readline()

for kpoint in kpoints:
    for x in range(192):
        actual_line = OUTCAR_file.readline()
        energy_value = actual_line.split()[1]

        kpoint[x] = energy_value

    for x in range(3):
        OUTCAR_file.readline()

OUTCAR_file.close()

band_matrix = np.zeros((192, 11))
y = 0

points_ordered = [1, 2, 3, 5, 7, 4, 1, 8, 11, 9, 3]

for kpoint in points_ordered:
    for x in range(192):
        band_matrix[x,y] = kpoints[kpoint-1][x]
    y = y + 1




kspace1 = np.linspace(0, 1, 11)
interv1 = np.linspace(0, 1, 110)
k_space_interp = np.linspace(0, 1, 110)

band_matrix_interp = np.zeros((192, 110))

for x in range(192):
    band_matrix_interp[x,:] = interp1d(kspace1, band_matrix[x,:],
                                  kind='quadratic', fill_value='extrapolate')(interv1)


x_labels = ['$\\Gamma$', 'X', 'M', '$\\Gamma$', 'R', 'X']
x_ticks = [0, 0.2, 0.4, 0.6, 0.8, 1]
for x in range(5):
    axs[3,0].axvline(x_ticks[x], color='black', linestyle='--', linewidth=.8)
axs[3,0].set_xticks(ticks=x_ticks, labels=x_labels)

for x in range(192):
    if x <= 91:
        color_type = 'coral'
    else:
        color_type = 'slateblue'
    axs[3,0].plot(interv1, band_matrix_interp[x,:]-np.max(band_matrix[0:92,:]), color=color_type)

for x in range(92):
    for y in range(11):
        if band_matrix[x,y] == np.max(band_matrix[0:92,:]):
            ind_max = y
axs[3,0].plot(kspace1[ind_max], np.max(band_matrix[0:92,:])-np.max(band_matrix[0:92,:]), 'o', markersize=8, color='darkred')

for x in range(92,192):
    for y in range(11):
        if band_matrix[x,y] == np.min(band_matrix[92:-1,:]):
            ind_min = y
axs[3,0].plot(kspace1[ind_min], np.min(band_matrix[92:-1,:])-np.max(band_matrix[0:92,:]), 'o', markersize=8, color='navy')


path = 'data/' + material + '/TDOS.dat'
total_dos = np.loadtxt(path)
path = 'data/' + material + '/PDOS_Ag.dat'
Ag_dos = np.loadtxt(path)
path = 'data/' + material + '/PDOS_S.dat'
S_dos = np.loadtxt(path)
path = 'data/' + material + '/PDOS_I.dat'
Br_dos = np.loadtxt(path)

axs[3,1].plot(total_dos[:,-1], total_dos[:,0], color='black', linewidth=.9)
axs[3,1].fill_between(total_dos[:,-1], total_dos[:,0], color='grey', edgecolor='None', alpha=0.5)
axs[3,1].plot(Ag_dos[:,-1], Ag_dos[:,0], color='cornflowerblue', linewidth=.9)
axs[3,1].fill_between(Ag_dos[:,-1], Ag_dos[:,0], color='cornflowerblue', edgecolor='None', alpha=0.5)
axs[3,1].plot(S_dos[:,-1], S_dos[:,0], color='lime', linewidth=.9)
axs[3,1].fill_between(S_dos[:,-1], S_dos[:,0], color='lime', edgecolor='None', alpha=0.5)
axs[3,1].plot(Br_dos[:,-1], Br_dos[:,0], color='violet', linewidth=.9)
axs[3,1].fill_between(Br_dos[:,-1], Br_dos[:,0], color='violet', edgecolor='None', alpha=0.5)
axs[3,1].set_xlim(0,50)
axs[3,1].set_ylim(-5,7.5)
axs[3,1].set_yticklabels([])
major_locator = MultipleLocator(20) 
minor_locator = MultipleLocator(4)  

axs[3,1].xaxis.set_major_locator(major_locator)
axs[3,1].xaxis.set_minor_locator(minor_locator)




plt.tight_layout()
plt.savefig('bands2.pdf')