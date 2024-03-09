import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

dist = [0, 0.2, 0.4, 0.6]

Ag1_S = []
Ag2_S = []
Ag3_S = []

Pb_O1 = []
Pb_O2 = []
Pb_O3 = []

Ti_O1 = []
Ti_O2 = []
Ti_O3 = []

length_file = open('lengths_file.txt', 'r')
length_file.readline()
actual_line = length_file.readline()
Ag1_S.append(float(actual_line.split()[1]))
Ag1_S.append(float(actual_line.split()[2]))
Ag1_S.append(float(actual_line.split()[3]))
Ag1_S.append(float(actual_line.split()[4]))

actual_line = length_file.readline()
Ag2_S.append(float(actual_line.split()[1]))
Ag2_S.append(float(actual_line.split()[2]))
Ag2_S.append(float(actual_line.split()[3]))
Ag2_S.append(float(actual_line.split()[4]))

actual_line = length_file.readline()
Ag3_S.append(float(actual_line.split()[1]))
Ag3_S.append(float(actual_line.split()[2]))
Ag3_S.append(float(actual_line.split()[3]))
Ag3_S.append(float(actual_line.split()[4]))

actual_line = length_file.readline()
Pb_O1.append(float(actual_line.split()[1]))
Pb_O1.append(float(actual_line.split()[2]))
Pb_O1.append(float(actual_line.split()[3]))
Pb_O1.append(float(actual_line.split()[4]))

actual_line = length_file.readline()
Pb_O2.append(float(actual_line.split()[1]))
Pb_O2.append(float(actual_line.split()[2]))
Pb_O2.append(float(actual_line.split()[3]))
Pb_O2.append(float(actual_line.split()[4]))

actual_line = length_file.readline()
Pb_O3.append(float(actual_line.split()[1]))
Pb_O3.append(float(actual_line.split()[2]))
Pb_O3.append(float(actual_line.split()[3]))
Pb_O3.append(float(actual_line.split()[4]))

actual_line = length_file.readline()
Ti_O1.append(float(actual_line.split()[1]))
Ti_O1.append(float(actual_line.split()[2]))
Ti_O1.append(float(actual_line.split()[3]))
Ti_O1.append(float(actual_line.split()[4]))

actual_line = length_file.readline()
Ti_O2.append(float(actual_line.split()[1]))
Ti_O2.append(float(actual_line.split()[2]))
Ti_O2.append(float(actual_line.split()[3]))
Ti_O2.append(float(actual_line.split()[4]))

actual_line = length_file.readline()
Ti_O3.append(float(actual_line.split()[1]))
Ti_O3.append(float(actual_line.split()[2]))
Ti_O3.append(float(actual_line.split()[3]))
Ti_O3.append(float(actual_line.split()[4]))

length_file.close()

fig, ax = plt.subplots(figsize=(4,3))
ax.set_ylim(1,3)
ax.set_xlabel('u (Å)')
ax.set_ylabel('Bond length (Å)')
ax.plot(dist, Ag1_S, marker='o', linestyle='-', color='lightcoral', alpha=0.6, label='Ag1-S')
ax.plot(dist, Ag2_S, marker='o', linestyle='-', color='darkred', alpha=0.6, label='Ag2-S')
ax.plot(dist, Ag3_S, marker='o', linestyle='-', color='chocolate', alpha=0.6, label='Ag3-S')
ax.plot(dist, Pb_O1, marker='o', linestyle='-', color='mediumspringgreen', alpha=0.6, label='Pb-O1')
ax.plot(dist, Pb_O2, marker='o', linestyle='-', color='paleturquoise', alpha=0.6, label='Pb-O2')
ax.plot(dist, Pb_O3, marker='o', linestyle='-', color='cyan', alpha=0.6, label='Pb-O3')
ax.plot(dist, Ti_O1, marker='o', linestyle='-', color='violet', alpha=0.6, label='Ti-O1')
ax.plot(dist, Ti_O2, marker='o', linestyle='-', color='orchid', alpha=0.6, label='Ti-O2')
ax.plot(dist, Ti_O3, marker='o', linestyle='-', color='hotpink', alpha=0.6, label='Ti-O3')

ax.legend(fontsize=5)

major_locator = MultipleLocator(0.1)  
minor_locator = MultipleLocator(0.02) 
ax.xaxis.set_major_locator(major_locator)
ax.xaxis.set_minor_locator(minor_locator)

major_locator = MultipleLocator(.5)  
minor_locator = MultipleLocator(0.1) 
ax.yaxis.set_major_locator(major_locator)
ax.yaxis.set_minor_locator(minor_locator)

plt.tight_layout()
plt.savefig('bond-length.pdf')