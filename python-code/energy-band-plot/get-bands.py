import numpy as np 
import matplotlib.pyplot as plt 

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

path = 'step2-test/OUTCAR'

OUTCAR_file = open(path, "r")
for x in range(1767):
    OUTCAR_file.readline()

for kpoint in kpoints:
    for x in range(192):
        actual_line = OUTCAR_file.readline()
        energy_value = actual_line.split()[1]

        kpoint[x] = energy_value

    for x in range(3):
        OUTCAR_file.readline()

OUTCAR_file.close()

band_matrix = np.zeros((192, 20))
y = 0

points_ordered = [1, 2, 3, 4, 7, 9, 10, 8, 5, 1, 11, 17, 20, 18, 13, 4, 20, 19, 16, 10]

for kpoint in points_ordered:
    for x in range(192):
        band_matrix[x,y] = kpoints[kpoint-1][x]
    y = y + 1


from scipy.interpolate import interp1d

kspace1 = np.linspace(0, 1, 16)
kspace2 = np.linspace(1, 1.2, 4)

interv1 = np.linspace(0, 1, 80)
interv2 = np.linspace(1, 1.2, 16)
k_space_interp = np.linspace(0, 1.2, 96)

band_matrix_interp = np.zeros((192, 96))

for x in range(192):
    band_matrix_interp[x,0:80] = interp1d(kspace1, band_matrix[x,0:16],
                                  kind='quadratic', fill_value='extrapolate')(interv1)
    band_matrix_interp[x,80:96] = interp1d(kspace2, band_matrix[x,16:20],
                                  kind='quadratic', fill_value='extrapolate')(interv2)

plt.figure()
plt.xlim(0, 1.2)
plt.ylim(-7.5, 11.5)
plt.ylabel('Energy (eV)')

x_labels = ['$\\Gamma$', 'X', 'M', '$\\Gamma$', 'R', 'X|R', 'M']
x_ticks = [0, 0.2, 0.4, 0.6, 0.8, 1, 1.2]
for x in range(5):
    plt.axvline(x_ticks[x], color='black', linestyle='--', linewidth=.8)
plt.xticks(ticks=x_ticks, labels=x_labels)

for x in range(192):
    if x <= 22:
        color_type = 'coral'
    else:
        color_type = 'slateblue'
    plt.plot(interv1, band_matrix_interp[x,0:80]-np.max(band_matrix[0:22,:]), color=color_type)
    plt.plot(interv2, band_matrix_interp[x,80:96]-np.max(band_matrix[0:22,:]), color=color_type)

for x in range(22):
    for y in range(16):
        if band_matrix[x,y] == np.max(band_matrix[0:22,:]):
            ind_max = y
plt.plot(kspace1[ind_max], np.max(band_matrix[0:22,:])-np.max(band_matrix[0:22,:]), 'o', markersize=10, color='darkred')

for x in range(23,192):
    for y in range(16):
        if band_matrix[x,y] == np.min(band_matrix[23:-1,:]):
            ind_min = y
plt.plot(kspace1[ind_min], np.min(band_matrix[23:-1,:])-np.max(band_matrix[0:22,:]), 'o', markersize=10, color='navy')


plt.axvline(1, color='black', linewidth=1.2)

plt.savefig('energy-band.pdf')

print('The maximum value is: ', np.max(band_matrix[0:22,:]))
print('The minimum value is: ', np.min(band_matrix[23:-1,:]))
print('Thus, the band gap is: ', np.min(band_matrix[23:-1,:]) - np.max(band_matrix[0:22,:]))