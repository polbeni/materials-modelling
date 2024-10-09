# Pol Benítez Colominas, October 2024
# University of Cambridge and Universitat Politècnica de Catalunya

# Computes the Fröhlich correction to band gap following the description in https://doi.org/10.1103/PhysRevB.102.045126

import numpy as np
import matplotlib.pyplot as plt


def polaron_ct(dielectric_ct, static_ct, effective_mass, freq_LO):
    """
    Determines the dimensionless polaron constant

    Inputs:
        dielectric_ct -> value of the high-frequency permittivity
        static_ct -> value of the static relative permittivity
        effective_mass -> value of the effective mass of the carrier
        freq_LO -> value of the LO freq after LOTO-splitting
    """

    alpha = ((e_charge**2) / (4 * np.pi * epsilon_0 * reduced_planck)) * ((1 / dielectric_ct) - (1 / static_ct)) * (effective_mass / (2 * reduced_planck * freq_LO))**0.5

    return alpha

def q_LO(effective_mass, freq_LO, freq):
    """
    Determines q_LO term

    Inputs:
        effective_mass -> value of the effective mass of the carrier
        freq_LO -> value of the LO freq after LOTO-splitting
        freq -> value of the freq
    """

    value = (2 * effective_mass * (freq_LO + freq) / reduced_planck)**0.5

    return value

def q_F(lattice_parameter, supercell_size):
    """
    Determines the q_F term

    Inputs:
        lattice_parameter -> lattice parameter of the unit cell (cubic cell assumed)
        supercell_size -> size of the supercell
    """

    volume_reciprocal = ((2 * np.pi) / (lattice_parameter * supercell_size))**3
    radius = ((3 * volume_reciprocal) / (4 * np.pi))**(1 / 3)
    
    return radius

def bose_eisntein_occupation(freq_LO, temp):
    """
    Determines the Bose-Einstein occupation factor of the LO mode

    Inputs:
        freq_LO -> value of the LO freq after LOTO-splitting
        temp -> temperature
    """

    value = 1 / (np.exp((reduced_planck * freq_LO) / (boltzmann_ct * temp)) - 1)

    return value

# define the basic constants from our DFT calculations
freq = 2 * np.pi * 3e12 #8e12 # THz 
freq_LO = 2 * np.pi * 4e12 #8.5e12 # THz
dielectric_ct = 7.182157 # adimensional
static_ct = 120.872611 + 7.182157 # adimensional
effective_mass_e = 0.842 * 9.11e-31 # electron effective mass, kg
effective_mass_h = 2.542 * 9.11e-31 # electron effective mass, kg
supercell_size = 4
lattice_value = 4.79338e-10 # in m 

# define the other constants we need
e_charge = 1.602e-19
epsilon_0 = 8.854e-12
reduced_planck = 1.054571817e-34
boltzmann_ct = 1.380649e-23

# let's compute everything
polaron_ct_value_e = polaron_ct(dielectric_ct, static_ct, effective_mass_e, freq_LO)
polaron_ct_value_h = polaron_ct(dielectric_ct, static_ct, effective_mass_h, freq_LO)
print(polaron_ct_value_e, polaron_ct_value_h)
q_F_value = q_F(lattice_value, supercell_size)
q_LO_value_e = q_LO(effective_mass_e, freq_LO, freq)
q_LO_value_h = q_LO(effective_mass_h, freq_LO, freq)

temperature = np.linspace(0.1, 1000, 10000)
bg_corr = np.zeros(10000)

it_value = 0
for t in temperature:
    be = bose_eisntein_occupation(freq_LO, t)

    vb_corr = (2 / np.pi) * polaron_ct_value_h * reduced_planck * freq_LO * (1 / np.tan(q_F_value / q_LO_value_h)) * (2 * be + 1)
    cb_corr = (2 / np.pi) * polaron_ct_value_e * reduced_planck * freq_LO * (1 / np.tan(q_F_value / q_LO_value_e)) * (2 * be + 1)

    bg_corr[it_value] = (cb_corr - vb_corr) * 6.242e21 # change to meV

    it_value = it_value + 1

# plot the results
fig, axs = plt.subplots(figsize=(4, 3))

axs.set_xlabel('T (K)')
axs.set_ylabel('$E_{g}$ Fröhlich correction (meV)')
axs.set_xlim(0, 1000)
axs.plot(temperature, bg_corr, color='black')

plt.tight_layout()
plt.savefig('frohlich_correction.pdf')