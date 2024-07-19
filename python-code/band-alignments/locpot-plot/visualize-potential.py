import matplotlib.pyplot as plt

distance = [] # in angstrom
potential = [] # planar average potential in eV

with open('PLANAR_AVERAGE.dat', 'r') as file:
    next(file) # jump the first line

    for line in file:
        values = line.split()

        distance.append(float(values[0]))
        potential.append(float(values[1]))

plt.figure()
plt.xlabel('Distance ($\AA$)')
plt.ylabel('Potential (eV)')
plt.plot(distance, potential)
plt.show()