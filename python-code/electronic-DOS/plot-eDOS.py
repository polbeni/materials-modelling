import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from scipy.interpolate import interp1d

def mean_val(min_x, max_x, num_points, array_curves):
    """
    This function computes the mean value of a curve and return the
    mean curve and the statistical error for each point

    Inputs: 
        min_x: minimum value of the interval for the x axis
        max_x: maximum value of the interval for the x axis
        num_points: number of points that we want for the interval
        array_curves_x: array that contains the x values for each curve
        array_curves_y: array that contains the y values for each curve
    Outputs:
        mean_curve: the computed mean curve
        error_curve: the error of the curve
        interval_curve: the x values of the curve
    """

    interval_curve = np.linspace(min_x, max_x, num_points)
    mean_curve = np.zeros(num_points)
    error_curve = np.zeros(num_points)

    for curve in range(len(array_curves)):
        curve_interpol = interp1d(array_curves[curve][:,0], array_curves[curve][:,-1],
                                  kind='linear', fill_value='extrapolate')(interval_curve)
        
        for element in range(num_points):
            mean_curve[element] = mean_curve[element] + curve_interpol[element]
    
    mean_curve[:] = mean_curve[:]/len(array_curves)

    partial_sum = np.zeros(num_points)
    for curve in range(len(array_curves)):
        curve_interpol = interp1d(array_curves[curve][:,0], array_curves[curve][:,-1],
                                  kind='linear', fill_value='extrapolate')(interval_curve)
        
        partial_sum[:] = partial_sum[:] + (curve_interpol[:] - mean_curve[:])**2
    
    error_curve[:] = ((1/(len(array_curves) - 1))*partial_sum[:])**.5


    return mean_curve, error_curve, interval_curve




t_0_total = np.loadtxt('data/Ag3SBr/T0/TDOS_SOC.dat')
t_0_Ag = np.loadtxt('data/Ag3SBr/T0/PDOS_Ag_SOC.dat')
t_0_S = np.loadtxt('data/Ag3SBr/T0/PDOS_S_SOC.dat')
t_0_Br = np.loadtxt('data/Ag3SBr/T0/PDOS_Br_SOC.dat')


arr_200_total = [None]*9
arr_200_Ag = [None]*9
arr_200_S = [None]*9
arr_200_Br = [None]*9

Ag3SBr = ['simulation-001', 'simulation-003', 'simulation-004', 'simulation-005', 
          'simulation-006', 'simulation-007', 'simulation-008', 'simulation-009', 
          'simulation-010']

iteration = 0
for simulation in Ag3SBr:
    path = 'data/Ag3SBr/T200/' + simulation + '/TDOS_SOC.dat'
    arr_200_total[iteration] = np.loadtxt(path)

    path = 'data/Ag3SBr/T200/' + simulation + '/PDOS_Ag_SOC.dat'
    arr_200_Ag[iteration] = np.loadtxt(path)

    path = 'data/Ag3SBr/T200/' + simulation + '/PDOS_S_SOC.dat'
    arr_200_S[iteration] = np.loadtxt(path)

    path = 'data/Ag3SBr/T200/' + simulation + '/PDOS_Br_SOC.dat'
    arr_200_Br[iteration] = np.loadtxt(path)

    iteration = iteration + 1

t_200_total, error_200_total, interval_x = mean_val(-8, 7, 1000, arr_200_total)
t_200_Ag, error_200_Ag, interval_x = mean_val(-8, 7, 1000, arr_200_Ag)
t_200_S, error_200_S, interval_x = mean_val(-8, 7, 1000, arr_200_S)
t_200_Br, error_200_Br, interval_x = mean_val(-8, 7, 1000, arr_200_Br)



plt.figure()
fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, figsize=(5,8))
fig.suptitle('Ag$_3$SBr $Pm\overline{3}m$')

#ax1.set_title('T=0 K')
#ax1.set_xlabel('E-E$_{F}$ (eV)')
ax1.set_ylabel('PDOS')
ax1.set_xlim(-4,4)
ax1.set_ylim(0,10)
ax1.plot(t_0_total[:,0], t_0_total[:,-1], linewidth=1, color='black', label='Total')
ax1.plot(t_0_Ag[:,0], t_0_Ag[:,-1], linewidth=1, color='royalblue', label='Ag')
ax1.plot(t_0_S[:,0], t_0_S[:,-1], linewidth=1, color='salmon', label='S')
ax1.plot(t_0_Br[:,0], t_0_Br[:,-1], linewidth=1, color='lightgreen', label='Br')
ax1.legend()

#ax2.set_title('T=200 K')
#ax2.set_xlabel('E-E$_{F}$ (eV)')
ax2.set_ylabel('PDOS')
ax2.set_xlim(-4,4)
ax2.set_ylim(0,10)
ax2.plot(interval_x[:], t_200_total[:]/8, linewidth=1, color='black', label='Total')
ax2.fill_between(interval_x[:], (t_200_total[:]/8 - error_200_total[:]/8),
                 (t_200_total[:]/8 + error_200_total[:]/8), color='black', alpha=0.6)
ax2.plot(interval_x[:], t_200_Ag[:]/8, linewidth=1, color='royalblue', label='Ag')
ax2.fill_between(interval_x[:], (t_200_Ag[:]/8 - error_200_Ag[:]/8),
                 (t_200_Ag[:]/8 + error_200_Ag[:]/8), color='royalblue', alpha=0.6)
ax2.plot(interval_x[:], t_200_S[:]/8, linewidth=1, color='salmon', label='S')
ax2.fill_between(interval_x[:], (t_200_S[:]/8 - error_200_S[:]/8),
                 (t_200_S[:]/8 + error_200_S[:]/8), color='salmon', alpha=0.6)
ax2.plot(interval_x[:], t_200_Br[:]/8, linewidth=1, color='lightgreen', label='Br')
ax2.fill_between(interval_x[:], (t_200_Br[:]/8 - error_200_Br[:]/8),
                 (t_200_Br[:]/8 + error_200_Br[:]/8), color='lightgreen', alpha=0.6)
#ax2.legend()


arr_400_total = [None]*9
arr_400_Ag = [None]*9
arr_400_S = [None]*9
arr_400_Br = [None]*9

Ag3SBr = ['simulation-001', 'simulation-002', 'simulation-003', 'simulation-004', 
          'simulation-005', 'simulation-006', 'simulation-007', 'simulation-009', 
          'simulation-010']

iteration = 0
for simulation in Ag3SBr:
    path = 'data/Ag3SBr/T400/' + simulation + '/TDOS_SOC.dat'
    arr_400_total[iteration] = np.loadtxt(path)

    path = 'data/Ag3SBr/T400/' + simulation + '/PDOS_Ag_SOC.dat'
    arr_400_Ag[iteration] = np.loadtxt(path)

    path = 'data/Ag3SBr/T400/' + simulation + '/PDOS_S_SOC.dat'
    arr_400_S[iteration] = np.loadtxt(path)

    path = 'data/Ag3SBr/T400/' + simulation + '/PDOS_Br_SOC.dat'
    arr_400_Br[iteration] = np.loadtxt(path)

    iteration = iteration + 1

t_400_total, error_400_total, interval_x = mean_val(-8, 7, 1000, arr_400_total)
t_400_Ag, error_400_Ag, interval_x = mean_val(-8, 7, 1000, arr_400_Ag)
t_400_S, error_400_S, interval_x = mean_val(-8, 7, 1000, arr_400_S)
t_400_Br, error_400_Br, interval_x = mean_val(-8, 7, 1000, arr_400_Br)


#ax3.set_title('T=400 K')
#ax3.set_xlabel('E-E$_{F}$ (eV)')
ax3.set_ylabel('PDOS')
ax3.set_xlim(-4,4)
ax3.set_ylim(0,10)
ax3.plot(interval_x[:], t_400_total[:]/8, linewidth=1, color='black', label='Total')
ax3.fill_between(interval_x[:], (t_400_total[:]/8 - error_400_total[:]/8),
                 (t_400_total[:]/8 + error_400_total[:]/8), color='black', alpha=0.6)
ax3.plot(interval_x[:], t_400_Ag[:]/8, linewidth=1, color='royalblue', label='Ag')
ax3.fill_between(interval_x[:], (t_400_Ag[:]/8 - error_400_Ag[:]/8),
                 (t_400_Ag[:]/8 + error_400_Ag[:]/8), color='royalblue', alpha=0.6)
ax3.plot(interval_x[:], t_400_S[:]/8, linewidth=1, color='salmon', label='S')
ax3.fill_between(interval_x[:], (t_400_S[:]/8 - error_400_S[:]/8),
                 (t_400_S[:]/8 + error_400_S[:]/8), color='salmon', alpha=0.6)
ax3.plot(interval_x[:], t_400_Br[:]/8, linewidth=1, color='lightgreen', label='Br')
ax3.fill_between(interval_x[:], (t_400_Br[:]/8 - error_400_Br[:]/8),
                 (t_400_Br[:]/8 + error_400_Br[:]/8), color='lightgreen', alpha=0.6)
#ax3.legend()


arr_600_total = [None]*8
arr_600_Ag = [None]*8
arr_600_S = [None]*8
arr_600_Br = [None]*8

Ag3SBr = ['simulation-001', 'simulation-002', 'simulation-003', 'simulation-004', 
          'simulation-005', 'simulation-007', 'simulation-008', 
          'simulation-009']

iteration = 0
for simulation in Ag3SBr:
    path = 'data/Ag3SBr/T600/' + simulation + '/TDOS_SOC.dat'
    arr_600_total[iteration] = np.loadtxt(path)

    path = 'data/Ag3SBr/T600/' + simulation + '/PDOS_Ag_SOC.dat'
    arr_600_Ag[iteration] = np.loadtxt(path)

    path = 'data/Ag3SBr/T600/' + simulation + '/PDOS_S_SOC.dat'
    arr_600_S[iteration] = np.loadtxt(path)

    path = 'data/Ag3SBr/T600/' + simulation + '/PDOS_Br_SOC.dat'
    arr_600_Br[iteration] = np.loadtxt(path)

    iteration = iteration + 1

t_600_total, error_600_total, interval_x = mean_val(-8, 7, 1000, arr_600_total)
t_600_Ag, error_600_Ag, interval_x = mean_val(-8, 7, 1000, arr_600_Ag)
t_600_S, error_600_S, interval_x = mean_val(-8, 7, 1000, arr_600_S)
t_600_Br, error_600_Br, interval_x = mean_val(-8, 7, 1000, arr_600_Br)


#ax4.set_title('T=600 K')
ax4.set_xlabel('E-E$_{F}$ (eV)')
ax4.set_ylabel('PDOS')
ax4.set_xlim(-4,4)
ax4.set_ylim(0,10)
ax4.plot(interval_x[:], t_600_total[:]/8, linewidth=1, color='black', label='Total')
ax4.fill_between(interval_x[:], (t_600_total[:]/8 - error_600_total[:]/8),
                 (t_600_total[:]/8 + error_600_total[:]/8), color='black', alpha=0.6)
ax4.plot(interval_x[:], t_600_Ag[:]/8, linewidth=1, color='royalblue', label='Ag')
ax4.fill_between(interval_x[:], (t_600_Ag[:]/8 - error_600_Ag[:]/8),
                 (t_600_Ag[:]/8 + error_600_Ag[:]/8), color='royalblue', alpha=0.6)
ax4.plot(interval_x[:], t_600_S[:]/8, linewidth=1, color='salmon', label='S')
ax4.fill_between(interval_x[:], (t_600_S[:]/8 - error_600_S[:]/8),
                 (t_600_S[:]/8 + error_600_S[:]/8), color='salmon', alpha=0.6)
ax4.plot(interval_x[:], t_600_Br[:]/8, linewidth=1, color='lightgreen', label='Br')
ax4.fill_between(interval_x[:], (t_600_Br[:]/8 - error_600_Br[:]/8),
                 (t_600_Br[:]/8 + error_600_Br[:]/8), color='lightgreen', alpha=0.6)
#ax4.legend()

ax1.text(-.5, 8.5, 'T=0 K', fontsize=14)
ax2.text(-.75, 8.5, 'T=200 K', fontsize=14)
ax3.text(-.75, 8.5, 'T=400 K', fontsize=14)
ax4.text(-.75, 8.5, 'T=600 K', fontsize=14)

major_locator = MultipleLocator(1) 
minor_locator = MultipleLocator(.2) 
ax1.xaxis.set_major_locator(major_locator)
ax1.xaxis.set_minor_locator(minor_locator)
ax2.xaxis.set_major_locator(major_locator)
ax2.xaxis.set_minor_locator(minor_locator)
ax3.xaxis.set_major_locator(major_locator)
ax3.xaxis.set_minor_locator(minor_locator)
ax4.xaxis.set_major_locator(major_locator)
ax4.xaxis.set_minor_locator(minor_locator)

major_locator = MultipleLocator(2.5) 
minor_locator = MultipleLocator(.5) 
ax1.yaxis.set_major_locator(major_locator)
ax1.yaxis.set_minor_locator(minor_locator)
ax2.yaxis.set_major_locator(major_locator)
ax2.yaxis.set_minor_locator(minor_locator)
ax3.yaxis.set_major_locator(major_locator)
ax3.yaxis.set_minor_locator(minor_locator)
ax4.yaxis.set_major_locator(major_locator)
ax4.yaxis.set_minor_locator(minor_locator)

plt.tight_layout()
plt.savefig('eDOS-Ag3SBr.pdf')