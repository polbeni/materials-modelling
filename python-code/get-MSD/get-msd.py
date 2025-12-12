# Pol Benítez Colominas, December 2025
# The University of Tokyo and Universitat Politècnica de Catalunya

# Get the MSD from a MD simulation by using kinisi module

import matplotlib.pyplot as plt
from scipy import stats

from kinisi.analyze import DiffusionAnalyzer

def plot_msd(msd_data, color, output_name):
    """
    Plots the MSD results and shows the diffusion coefficient
    
    :param msd_data: kinsi object with the msd information
    :param color: color of the MSD curve and its error
    :param output_name: name of the output file
    """

    # Compute the diffusion coefficient
    n_points = len(msd.dt)
    start_idx = n_points // 2

    slope, _, _, _, _ = stats.linregress(msd.dt[start_idx:], msd.msd[start_idx:])

    D = slope / 6  # in Å^2/ps, we divide between 6 for the 3D diffusion

    # Plot the results of MSD and difussion coefficient
    _, axs = plt.subplots(figsize=(4, 3))

    axs.set_xlabel('$\Delta t$ (ps)')
    axs.set_ylabel('MSD (Å$^{2}$)')

    axs.plot(msd_data.dt, msd_data.msd, linewidth=2.5, color=color)

    axs.fill_between(msd_data.dt, msd_data.msd - msd_data.msd_std, msd_data.msd + msd_data.msd_std, alpha=0.3, color=color, edgecolor='none')

    axs.set_xlim(0, max(msd_data.dt))
    axs.set_ylim(0, max(msd_data.msd + msd_data.msd_std) + max(msd_data.msd + msd_data.msd_std)*0.1)

    textstr = f'D = {D:.3f} Å$^{2}$/ps'
    axs.text(0.05, 0.95, textstr, transform=axs.transAxes, fontsize=10,
            verticalalignment='top', horizontalalignment='left',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    plt.tight_layout()
    plt.savefig(output_name)


# Define your parameters
XDATCAR_path = 'XDATCAR'
params = {'specie': 'Au', 
          'time_step': 1.0, 
          'step_skip': 1}

# Obtain the MSD
msd = DiffusionAnalyzer.from_file(XDATCAR_path, parser_params=params)

# Save the results in a plot
plot_msd(msd, 'gold', 'Au_MSD.pdf')
