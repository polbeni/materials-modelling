# Pol Benítez Colominas, October 2024
# University of Cambridge and Universitat Politècnica de Catalunya

# Computes the Fröhlich correction to band gap following the description in https://doi.org/10.1103/PhysRevB.102.045126

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator


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
freq = [2 * np.pi * 0.4e12, 2 * np.pi * 1.82e12, 2 * np.pi * 7.88e12] # THz 
freq_LO = [2 * np.pi * 0.95e12, 2 * np.pi * 4.02e12, 2 * np.pi * 8.45e12] # THz
dielectric_ct = 7.182157 # adimensional
static_ct = 120.872611 + 7.182157 # adimensional
effective_mass_e = 0.842 * 9.11e-31 # electron effective mass, kg
effective_mass_h = 2.542 * 9.11e-31 # electron effective mass, kg
supercell_size = 4 # adimensional
lattice_value = 4.79338e-10 # in m 

# define the other constants we need
e_charge = 1.602e-19 # A·s
epsilon_0 = 8.854e-12 # A^2·kg^-1·m^-3·s^5
reduced_planck = 1.054571817e-34 # kg·m^2·s^-1
boltzmann_ct = 1.380649e-23 # kg^2·m^2·s^-2·K^-1

# compute the Fröhlich correction for each splitting
temperature = np.linspace(0.1, 1000, 10000)
bg_correction = []

for splitting in range(len(freq)):
    polaron_ct_value_e = polaron_ct(dielectric_ct, static_ct, effective_mass_e, freq_LO[splitting])
    polaron_ct_value_h = polaron_ct(dielectric_ct, static_ct, effective_mass_h, freq_LO[splitting])
    q_F_value = q_F(lattice_value, supercell_size)
    q_LO_value_e = q_LO(effective_mass_e, freq_LO[splitting], freq[splitting])
    q_LO_value_h = q_LO(effective_mass_h, freq_LO[splitting], freq[splitting])

    bg_corr = np.zeros(10000)

    it_value = 0
    for t in temperature:
        be = bose_eisntein_occupation(freq_LO[splitting], t)

        vb_corr = (2 / np.pi) * polaron_ct_value_h * reduced_planck * freq_LO[splitting] * (1 / np.tan(q_F_value / q_LO_value_h)) * (2 * be + 1)
        cb_corr = (2 / np.pi) * polaron_ct_value_e * reduced_planck * freq_LO[splitting] * (1 / np.tan(q_F_value / q_LO_value_e)) * (2 * be + 1)

        bg_corr[it_value] = (cb_corr - vb_corr) * 6.242e21 # change to meV

        it_value = it_value + 1

    bg_correction.append(bg_corr[:])

# correction with averaged frequencies
freq_mean = np.mean(freq)
freq_LO_mean = np.mean(freq_LO)

polaron_ct_value_e = polaron_ct(dielectric_ct, static_ct, effective_mass_e, freq_LO_mean)
polaron_ct_value_h = polaron_ct(dielectric_ct, static_ct, effective_mass_h, freq_LO_mean)
q_F_value = q_F(lattice_value, supercell_size)
q_LO_value_e = q_LO(effective_mass_e, freq_LO_mean, freq_mean)
q_LO_value_h = q_LO(effective_mass_h, freq_LO_mean, freq_mean)

bg_corr_mean = np.zeros(10000)

it_value = 0
for t in temperature:
    be = bose_eisntein_occupation(freq_LO_mean, t)

    vb_corr = (2 / np.pi) * polaron_ct_value_h * reduced_planck * freq_LO_mean * (1 / np.tan(q_F_value / q_LO_value_h)) * (2 * be + 1)
    cb_corr = (2 / np.pi) * polaron_ct_value_e * reduced_planck * freq_LO_mean * (1 / np.tan(q_F_value / q_LO_value_e)) * (2 * be + 1)

    bg_corr_mean[it_value] = (cb_corr - vb_corr) * 6.242e21 # change to meV

    it_value = it_value + 1


# plot the results
fig, axs = plt.subplots(figsize=(4, 2.6))

axs.set_xlabel('T (K)')
axs.set_ylabel('Fröhlich $E_{g}$ correction (meV)')
axs.set_xlim(0, 600)
axs.set_ylim(-300, 50)
axs.plot(temperature, bg_correction[0], color='black', label='0.40 THz $\\to$ 0.95 THz')
axs.plot(temperature, bg_correction[1], color='grey', label='1.82 THz $\\to$ 4.02 THz')
axs.plot(temperature, bg_correction[2], color='blue', label='7.88 THz $\\to$ 8.45 THz')

axs.plot(temperature, bg_corr_mean, color='red', label='$\overline{\\omega}$ $\\to$ $\overline{\\omega}_{LO}$')

axs.legend(frameon=False, fontsize=8)


major_locator = MultipleLocator(50) 
minor_locator = MultipleLocator(10)  

axs.yaxis.set_major_locator(major_locator)
axs.yaxis.set_minor_locator(minor_locator)

major_locator = MultipleLocator(100) 
minor_locator = MultipleLocator(20)  

axs.xaxis.set_major_locator(major_locator)
axs.xaxis.set_minor_locator(minor_locator)


plt.tight_layout()
plt.savefig('frohlich_correction-multiple-splitting.pdf')