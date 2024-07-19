# Pol Benítez Colominas, July 2024
# Universitat Politècnica de Catalunya

# Generates the band alignment from bulk and slab LOCPOT calculations
# extract information from LOCPOT using VASPKIT, and put PLANAR_AVERAGE.dat
# files in the corresponding path (use option 426)

import subprocess
import numpy as np

def electronic_bandGap(file_name):
    """
    This functions uses DOSCAR file generated in VASP simulations and returns the Fermi energy
    the band gap, and the energies of the band gap (respect the exchange-correlation functional
    used).

    Inputs:
        file_name: path of the DOSCAR file
    """
    
    file = open(file_name, "r")

    for x in range(6):
        actual_string = file.readline()
        if x == 5:
            fermiEnergy = float(actual_string.split()[3])

    file.close()

    file = open(file_name, "r")

    for x in range(6):
        file.readline()

    for x in file:
        actual_string = x

        if (float(actual_string.split()[0]) <= fermiEnergy+0.1) and (float(actual_string.split()[0]) >= fermiEnergy-0.1):
            density_bandGap = float(actual_string.split()[2])

            break

    file.close()

    file = open(file_name, "r")

    for x in range(6):
        file.readline()

    for x in file:
        actual_string = x

        if float(actual_string.split()[2]) == density_bandGap:
            minEnergy = float(actual_string.split()[0])

            break   

    for x in file:
        actual_string = x

        if float(actual_string.split()[2]) != density_bandGap:
            maxEnergy = float(actual_string.split()[0])

            break 
    bandGap = maxEnergy - minEnergy

    file.close()
    
    return fermiEnergy, minEnergy, maxEnergy, bandGap

def average_material_slab(potential):
    """
    Computes the average potential for slab ignoring the refered points to vacuum
    (it is done for the Pm-3m phase where we have a potential with two kinds of minimums
    in different systems some parameters should be changed)

    Inputs:
        potential: values of the potential for the slab structure
    """

    num_max_minimums = 4

    points = []

    num_max = 0
    iteration = 120
    val1 = 1000
    first_condition = False
    while first_condition == False:
        val2 = potential[iteration]
        if val2 > val1:
            first_condition = True
        else:
            val1 = val2
            iteration = iteration + 1

    up = True
    while num_max < num_max_minimums:
        val2 = potential[iteration]
        points.append(val2)

        if up == True:
            if val2 < val1:
                up = False
                num_max = num_max + 1
        else:
            if val2 > val1:
                up = True
                num_max = num_max + 1

        val1 = val2
        iteration = iteration + 1

    points = np.array(points)
    value = np.mean(points)

    return value

materials = ['Ag3SCl', 'Ag3SBr', 'Ag3SI', 'Ag3SeCl', 'Ag3SeBr', 'Ag3SeI']

for material in materials:
    path_results = 'data/' + material + '/'

    # potential bulk
    distance_bulk = [] # in angstrom
    potential_bulk = [] # planar average potential in eV

    with open(path_results + 'bulk/PLANAR_AVERAGE.dat', 'r') as file:
        next(file) # jump the first line

        for line in file:
            values = line.split()

            distance_bulk.append(float(values[0]))
            potential_bulk.append(float(values[1]))

    distance_bulk = np.array(distance_bulk)
    potential_bulk = np.array(potential_bulk)

    # potential slab
    distance_slab = [] # in angstrom
    potential_slab = [] # planar average potential in eV

    with open(path_results + 'slab/PLANAR_AVERAGE.dat', 'r') as file:
        next(file) # jump the first line

        for line in file:
            values = line.split()

            distance_slab.append(float(values[0]))
            potential_slab.append(float(values[1]))

    distance_slab = np.array(distance_slab)
    potential_slab = np.array(potential_slab)

    # define the indices of the vacuum line in slab (in our case we have vacuum-slab-vacuum,
    # then we need the indices of the vacuum at the left and the right of the slab)
    vacuum_indices_left = range(0, int(len(distance_slab)/10))
    vacuum_indices_right = range(len(distance_slab) - int(len(distance_slab)/10), len(distance_slab))

    # compute the vacuum level
    vacuum_level = np.mean([np.mean(potential_slab[vacuum_indices_left]), np.mean(potential_slab[vacuum_indices_right])])

    # compute the average potential of the bulk and slab
    bulk_average_potential = np.mean(potential_bulk)
    slab_average_potential = average_material_slab(potential_slab)

    # define the valence band (VB) and the conduction band (CB) level
    _, VB, CB, _ = electronic_bandGap(path_results + 'bulk/DOSCAR')

    # compute the absolute bands
    VB_absolute = VB + (slab_average_potential - bulk_average_potential ) - vacuum_level
    CB_absolute = CB + (slab_average_potential - bulk_average_potential ) - vacuum_level

    print(f'Band alignment results for {material}')
    print(f'VB: {VB_absolute} eV')
    print(f'CB: {CB_absolute} eV')
    print('')