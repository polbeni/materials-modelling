# Pol Benítez Colominas, June 2025
# Universitat Politècnica de Catalunya

# Computes the enthalpy given a E(V) list

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, root_scalar

def birch_murnaghan(V, E0, V0, B0, B0_prime):
    """
    Birch-Murnaghan EOS

    Inputs:
        V: volume (Angstrom**3)
        E0, V0, B0, B0_prime: Birch-Murnaghan parameters
    """

    E = E0 + ((9 * V0 * B0) / (16)) * ((((V0 / V)**(2/3) - 1)**3 * B0_prime) + (((V0/V)**(2/3) - 1)**2 * (6 - 4*(V0/V)**(2/3))))

    return E

def fit_EOS(volume_list, energy_list):
    """
    Fit the desired EOS (Equation Of State)

    Inputs:
        volume_list: list of volumes (Angstrom**3)
        energy_list: list of energies (eV)
    """

    volume_list = np.array(volume_list)
    energy_list = np.array(energy_list)

    # Give an initial guess for the parameters
    E0_guess = min(energy_list)
    V0_guess = volume_list[np.argmin(energy_list)]
    B0_guess = 0.5
    B0_prime_guess = 4.0

    # Fit the parameters
    popt, _ = curve_fit(birch_murnaghan, volume_list, energy_list, p0=[E0_guess, V0_guess, B0_guess, B0_prime_guess])

    return popt

def pressure_birch_murnaghan(V, E0, V0, B0, B0_prime):
    """
    Pressure determined from the Birch-Murnaghan EOS

    Inputs:
        V: volume (Angstrom**3)
        E0, V0, B0, B0_prime: Birch-Murnaghan parameters
    """

    pressure = ((3 * B0) / 2) * ((V0 / V)**(7/3) - (V0 / V)**(5/3)) * (1 + (3 / 4) * (B0_prime - 4) * ((V0 / V)**(2/3) - 1))

    return pressure

def determine_H_vs_p(E0, V0, B0, B0_prime, pressure_list):
    """
    Determines the enthalpy vs pressure relation, H(p), from the E(V) data and the Birch-Murnaghan EOS

    Inputs:
        E0, V0, B0, B0_prime: Birch-Murnaghan parameters
        pressure_list: list with the desired pressures to compute the enthalpy
    """

    # Define functions to determine E(V) and p(V) relations
    def E(V): return birch_murnaghan(V, E0, V0, B0, B0_prime)
    def p(V): return pressure_birch_murnaghan(V, E0, V0, B0, B0_prime)

    # Create list for enthalpy and volumes
    V_list = []
    H_list = []

    # Calculate the elements
    for p_target in pressure_list:
        # Root finding to solve p(V) = p_target
        sol = root_scalar(lambda V: p(V) - p_target, bracket=[0.5 * V0, 1.5 * V0], method='brentq')

        if sol.converged:
            V = sol.root
            H = E(V) + p_target * V

            V_list.append(V)
            H_list.append(H)
        else:
            V_list.append(np.nan)
            H_list.append(np.nan)

    return V_list, H_list


# E(V) data points
v_data = [130.5, 132, 133.5, 135, 136.2, 137.7, 139.4, 140.7, 141.8]
e_data = [-121.75, -121.80, -121.83, -121.84, -121.83, -121.79, -121.76, -121.71, -121.65]

# Fit the parameters
parameters = fit_EOS(v_data, e_data)
print(parameters)

# fit the Birch-Murnaghan 
volumes = np.linspace(130, 142, 1000)
energies = np.zeros(1000)
for x in range(len(volumes)):
    energies[x] = birch_murnaghan(volumes[x], parameters[0], parameters[1], parameters[2], parameters[3])

# Plot the E(V) and fitted Birch-Murnaghan
plt.figure()
plt.xlabel('V ($\AA^{3}$)')
plt.ylabel('E (eV)')
plt.plot(v_data, e_data, marker='o', linestyle='', label='points')
plt.plot(volumes, energies, label='fitted EOS')
plt.legend()
plt.tight_layout()
plt.savefig('E_vs_V.pdf')

# Determine the H(p) relation
min_pressure = 0 # in GPa
max_pressure = 100 # in GPa

pressure_array = np.linspace(min_pressure, max_pressure, 1000) # in GPa

min_pressure = min_pressure * 0.06241509 # in eV/A^3
max_pressure = max_pressure * 0.06241509 # in eV/A^3

pressure_array_evA = np.linspace(min_pressure, max_pressure, 1000) # in eV/A^3

_, enthalpy = determine_H_vs_p(parameters[0], parameters[1], parameters[2], parameters[3], pressure_array_evA)

# Plot the enthalpy vs pressure relation
plt.figure()
plt.xlabel('p (GPa)')
plt.ylabel('H (eV)')
plt.plot(pressure_array, enthalpy)
plt.tight_layout()
plt.savefig('H_vs_p.pdf')