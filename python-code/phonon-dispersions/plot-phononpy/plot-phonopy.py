# Pol Benítez Colominas, August 2025
# Universitat Politècnica de Catalunya

# Plot the phonon dispersions from the band.yaml output file obtained using Phonopy software


### Modules
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import yaml


### Functions
def phonons_from_phonopy(file_path):
    """
    Extract data from band.yaml file

    Inputs:
        file_path: the path and name of the file

    Outputs:
        num_atoms: the number of atoms in the unit cell
        nqpoints: number of points in the reciprocal space
        npaths: the number of different paths in the reciprocal space
        segments_nqpoints: the number of points in the reciprocal space for each path
        points_labels: the labels of the high symmetry points
        phonons: a matrix array with the information about the phonons, it has the following structure
                    q-point   | branch 1 | branch 2 | branch 3 | ...
                    ------------------------------------------------
                    q-point 1 |   freq   |   freq   |   freq   | ...
                    q-point 2 |   freq   |   freq   |   freq   | ...
                    q-point 3 |   freq   |   freq   |   freq   | ...
                    ...       |   ...    |   ...    |   ...    | ...
    """

    with open(file_path, "r") as file:
        data = yaml.safe_load(file)

    num_atoms = data["natom"]
    nqpoints = data["nqpoint"]
    npaths = data["npath"]
    segments_nqpoints = data["segment_nqpoint"]
    points_labels = data["labels"]

    phonons = np.zeros((nqpoints, num_atoms*3 + 1))
    for x in range(nqpoints):
        phonons[x,0] = data["phonon"][x]["distance"]
        for y in range(num_atoms*3):
            phonons[x,y+1] = data["phonon"][x]["band"][y]["frequency"]

    return num_atoms, nqpoints, npaths, segments_nqpoints, points_labels, phonons

def plot_phonon(paths, colors, labels, want_legend, title, want_title, y_lims, ticks_num, alpha_disp, width_disp, output_name, threshold_disp):
    """
    Plots the phonon dispersions from band.yaml. It can handle different phonon data,
    but the same k-path should be provided

    Inputs:
        paths: path or paths to the band.yaml files (array)
        colors: list with the colors of the dispersions for each file (array)
        labels: list with the labels for each file (array)
        want_legend: indicate if you want to plot the legend (boolean)
        title: title for the plot (str)
        want_title: indicate if you want to show a title (boolean)
        y_lims: min and max value of y to show in the plot (array)
        alpha_disp: desired alpha value for the dispersions (float)
        width_disp: desired line width for the dispersions (float)
        ticks_num: list with the number for every big and small ticks (array)
        output_name: name for the output file
        threshold_disp: threshold for non continuous bands !!!(not really sure why, but take a small threshold, such as 1e-13)
    """

    plt.figure()
    _, axs = plt.subplots(figsize=(4, 3))

    if want_title == True:
        axs.set_title(title) #, pad=20)

    axs.set_ylabel('$\omega$ (THz)')

    if y_lims[0] < 0:
        axs.axhspan(ymin=y_lims[0]*2, ymax=0, color='grey', alpha=0.5)

    for file_phon in range(len(paths)):
        file_path = paths[file_phon]

        num_atoms, nqpoints, npaths, segments_nqpoints, points_labels, phonons = phonons_from_phonopy(file_path)

        if file_phon == 0:
            vertical_lines = ['False']*(npaths-1)
            for x in range(npaths-1):
                if points_labels[x][1] != points_labels[x+1][0]:
                    vertical_lines[x] = True
                    
            x_labels = ['point']*(npaths+1)
            for x in range(npaths):
                if x == 0:
                    x_labels[0] = points_labels[0][0]
                if  x+1 < npaths:
                    if vertical_lines[x] == True:
                        x_labels[x+1] = points_labels[x][1] + '|' + points_labels[x+1][0]
                    else: 
                        x_labels[x+1] = points_labels[x][1]
                else:
                    x_labels[x+1] = points_labels[x][1]
            
            axs.set_xlim((phonons[0,0],phonons[-1,0]))

            axs.axhline(0, color='black', linestyle='--', linewidth=1)

            k_point_number = 0
            num_segment = 0
            distance = 0
            for x in segments_nqpoints:
                k_point_number = k_point_number + x
                if vertical_lines[num_segment] == True:
                    distance = phonons[k_point_number-1,0]
                    axs.axvline(distance, color='black', linewidth=1)
                num_segment = num_segment + 1
                if num_segment >= npaths-1:
                    break

            x_ticks = [0]*(npaths+1)
            k_point_number = 0
            element_ticks = 1
            for x in segments_nqpoints:
                k_point_number = k_point_number + x
                x_ticks[element_ticks] = phonons[k_point_number-1,0]
                element_ticks = element_ticks + 1

            axs.set_xticks(ticks=x_ticks, labels=x_labels)

        k_path = phonons[:, 0]
        disp_continuity = []
        k_path_continuity = []
        for x in range(num_atoms*3):
            current_disp_curve = []
            current_disp_kpath = []

            actual_disp_line = phonons[:, x+1]

            total_points = 0

            fragment = []
            fragment_k = []
            for num_points in range(segments_nqpoints[0]):
                fragment.append(actual_disp_line[total_points])
                fragment_k.append(k_path[total_points])

                total_points = total_points + 1

            for segment_intersection in range(npaths - 1):
                if abs(actual_disp_line[total_points] - actual_disp_line[total_points + 1]) > threshold_disp:
                    current_disp_curve.append(fragment)
                    current_disp_kpath.append(fragment_k)

                    fragment = []
                    fragment_k = []

                for num_points in range(segments_nqpoints[segment_intersection + 1]):
                    fragment.append(actual_disp_line[total_points])
                    fragment_k.append(k_path[total_points])

                    total_points = total_points + 1
            
            current_disp_curve.append(fragment)
            current_disp_kpath.append(fragment_k)
            
            disp_continuity.append(current_disp_curve)
            k_path_continuity.append(current_disp_kpath)

        initial_disp_curve_ploted = True
        for x in range(num_atoms*3):
            for fragment in range(len(disp_continuity[x])):
                if (x == 0) and (initial_disp_curve_ploted == True):
                    axs.plot(k_path_continuity[x][fragment], disp_continuity[x][fragment], linestyle='solid', color=colors[file_phon], alpha=alpha_disp, linewidth=width_disp, label=labels[file_phon])

                    initial_disp_curve_ploted = False
                else:
                    axs.plot(k_path_continuity[x][fragment], disp_continuity[x][fragment], linestyle='solid', color=colors[file_phon], alpha=alpha_disp, linewidth=width_disp)
    
    major_locator = MultipleLocator(ticks_num[0]) 
    minor_locator = MultipleLocator(ticks_num[1])  
    axs.yaxis.set_major_locator(major_locator)
    axs.yaxis.set_minor_locator(minor_locator)

    axs.set_ylim(y_lims[0], y_lims[1])

    if want_legend == True:
        axs.legend(frameon=False, bbox_to_anchor=(0.5, 0.97), ncol=2, loc='lower center')

    plt.tight_layout()
    plt.savefig(output_name)


# Test the code
plot_phonon(['cubic/harmonic/band.yaml', 'cubic/anharmonic/band.yaml'], ['black', 'red'], ['harmonic', 'anharmonic'], False,
            '$Fm\overline{3}m$', True, [-6, 25], [5, 1], 0.6, 1, 'phonon-cubic.pdf', 1e-13)
