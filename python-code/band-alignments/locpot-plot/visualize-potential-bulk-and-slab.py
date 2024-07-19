import matplotlib.pyplot as plt

distance_bulk = [] # in angstrom
potential_bulk = [] # planar average potential in eV

with open('bulk/PLANAR_AVERAGE.dat', 'r') as file:
    next(file) # jump the first line

    for line in file:
        values = line.split()

        distance_bulk.append(float(values[0]))
        potential_bulk.append(float(values[1]))

distance_slab = [] # in angstrom
potential_slab = [] # planar average potential in eV

with open('slab/PLANAR_AVERAGE.dat', 'r') as file:
    next(file) # jump the first line

    for line in file:
        values = line.split()

        distance_slab.append(float(values[0]))
        potential_slab.append(float(values[1]))

plt.figure()
plt.xlabel('Distance ($\AA$)')
plt.ylabel('Potential (eV)')
plt.plot(distance_bulk, potential_bulk, label='bulk')
plt.plot(distance_slab, potential_slab, label='slab')
plt.legend()
plt.show()