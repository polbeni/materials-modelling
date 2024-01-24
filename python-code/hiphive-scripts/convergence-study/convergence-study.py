# Pol Benítez Colominas, January 2024
# Universitat Politècnica de Catalunya

# Code to study the convergence after using parameters.py script
# For now the only parameters studied are the R2 of train and test set as well as the 
# loss function metric of some high symmetry points in the reciprocal-space, after the 
# harmonic phonon frequencies determination

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

def loss_function(exact_freq, hiphive_freq):
    """
    Computes the loss function 

    Inputs:
        exact_freq: exact phonon frequencies at a given high symmetry point computed with Phonopy
        hiphive_freq: hiPhive determined phonon frequencies at a given high symmetry point
    """
    q_loss = 0

    for x in range(len(exact_freq)):
        q_loss = q_loss + (exact_freq[x] - hiphive_freq[x])**2

    return q_loss

### R2 convergence study

# create a file where storage the convergence results
results_r2 = open('convergence_r2.txt', 'w')
results_r2.write('cutoffs   #structures   R2_train   R2_test \n')

# define all the cutoffs of interest and number of structures
diff_cutoffs = [[4], [4.5], [5]]
diff_num_structures = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500]

# define an array with all the directories names
dir_results = ['1_050', '1_100', '1_150', '1_200', '1_250', '1_300', '1_350', '1_400', '1_450', '1_500',
               '2_050', '2_100', '2_150', '2_200', '2_250', '2_300', '2_350', '2_400', '2_450', '2_500',
               '3_050', '3_100', '3_150', '3_200', '3_250', '3_300', '3_350', '3_400', '3_450', '3_500']

num_iteration = 0

for cuts in diff_cutoffs:
    for struc in diff_num_structures:
        optm_file = open(dir_results[num_iteration] + '/log_files/4-optimizer', 'r')
        for x in range(11):
            optm_file.readline()

        actual_line = optm_file.readline()
        r2_train = actual_line.split()[2]
        actual_line = optm_file.readline()
        r2_test = actual_line.split()[2]

        results_r2.write(f'{cuts}   {struc}   {r2_train}   {r2_test} \n')

        optm_file.close()

        num_iteration = num_iteration + 1

results_r2.close()


### frequencies convergence study



### plot the results

fig, axs = plt.subplots(2, 1, figsize=(4,4))

axs[0].set_xlabel('# structures')
axs[0].set_ylabel('$R^2$ train')
axs[1].set_xlabel('# structures')
axs[1].set_ylabel('$R^2$ test')

num_structures = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500]

train_1 = []
train_2 = []
train_3 = []

test_1 = []
test_2 = []
test_3 = []

results_r2 = open('convergence_r2.txt', 'r')
results_r2.readline()

for x in range(len(num_structures)):
    actual_line = results_r2.readline()
    train_1.append(float(actual_line.split()[2]))
    test_1.append(float(actual_line.split()[3]))
    
for x in range(len(num_structures)):
    actual_line = results_r2.readline()
    train_2.append(float(actual_line.split()[2]))
    test_2.append(float(actual_line.split()[3]))

for x in range(len(num_structures)):
    actual_line = results_r2.readline()
    train_3.append(float(actual_line.split()[2]))
    test_3.append(float(actual_line.split()[3]))

results_r2.close()

axs[0].plot(num_structures, train_1, marker='o', linestyle='--', linewidth=0.5, 
            color='salmon', alpha=0.6, label='Cutoff: 4.0 Å')
axs[0].plot(num_structures, train_2, marker='o', linestyle='--', linewidth=0.5, 
            color='deepskyblue', alpha=0.6, label='Cutoff: 4.5 Å')
axs[0].plot(num_structures, train_3, marker='o', linestyle='--', linewidth=0.5, 
            color='lightgreen', alpha=0.6, label='Cutoff: 5.0 Å')

axs[1].plot(num_structures, test_1, marker='o', linestyle='--', linewidth=0.5, 
            color='salmon', alpha=0.6, label='Cutoff: 4.0 Å')
axs[1].plot(num_structures, test_2, marker='o', linestyle='--', linewidth=0.5, 
            color='deepskyblue', alpha=0.6, label='Cutoff: 4.5 Å')
axs[1].plot(num_structures, test_3, marker='o', linestyle='--', linewidth=0.5, 
            color='lightgreen', alpha=0.6, label='Cutoff: 5.0 Å')

axs[0].legend(bbox_to_anchor=(1.05, 1.0), loc='upper left', fontsize=7)

major_locator = MultipleLocator(0.015) 
minor_locator = MultipleLocator(0.003) 
axs[0].yaxis.set_major_locator(major_locator)
axs[0].yaxis.set_minor_locator(minor_locator)
axs[1].yaxis.set_major_locator(major_locator)
axs[1].yaxis.set_minor_locator(minor_locator)

major_locator = MultipleLocator(100) 
minor_locator = MultipleLocator(20) 
axs[0].xaxis.set_major_locator(major_locator)
axs[0].xaxis.set_minor_locator(minor_locator)
axs[1].xaxis.set_major_locator(major_locator)
axs[1].xaxis.set_minor_locator(minor_locator)

axs[0].set_ylim(0.92, 0.96)
axs[1].set_ylim(0.92, 0.96)

plt.tight_layout()
plt.savefig('r2-convergence.pdf')