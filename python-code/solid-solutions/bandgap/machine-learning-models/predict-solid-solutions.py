# Pol Benítez Colominas, July 2024
# Universitat Politècnica de Catalunya

# Script to use the machine learning trained model to predict band gaps for new solid solution compounds

import joblib

import numpy as np
import matplotlib.pyplot as plt

# Load the model
loaded_model = joblib.load('trained_model.joblib')

# Create the concentrations we want to predict
num_points = 1000
S_conc = np.linspace(0, 1, num_points)
Br_conc = np.linspace(0, 1, num_points)

solid_solutions = []
for x in range(len(S_conc)):
    for y in range(len(Br_conc)):
        row = [S_conc[x], Br_conc[y]]
        solid_solutions.append(row)

solid_solutions = np.array(solid_solutions)

# Predict the band gaps
y_pred = loaded_model.predict(solid_solutions)

# Generate the band gap matrix
matrix_bg = []
iteration = 0
for x in range(num_points):
    row = []
    for y in range(num_points):
        row.append(y_pred[(num_points**2 -1) - iteration])

        iteration = iteration + 1
    
    n_row = []
    for y in range(num_points):
        n_row.append(row[(num_points - 1) - y])
    
    matrix_bg.append(n_row)

# Plot the heat map for band gaps
fig, ax = plt.subplots(figsize=(6,4))

cax = ax.imshow(matrix_bg, cmap='magma', aspect='auto')
cbar = fig.colorbar(cax, ax=ax)
cbar.set_label('$E_{g}$ (eV)', fontsize=12) 

ax.set_title('Ag$_3$S$_{x}$Se$_{1-x}$Br$_{y}$I$_{1-y}$\n Machine Learning prediction')
ax.set_xlabel('x', fontsize=12)
ax.set_ylabel('y', fontsize=12)


x_ticks = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
for x in range(len(x_ticks)):
    x_ticks[x] = x_ticks[x]*(num_points/10)
x_labels = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 10]
ax.set_xticks(ticks=x_ticks, labels=x_labels)

y_ticks = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
for x in range(len(x_ticks)):
    y_ticks[x] = y_ticks[x]*(num_points/10)
y_labels = [1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0]
ax.set_yticks(ticks=y_ticks, labels=y_labels)


plt.tight_layout()
plt.savefig('solid-solutions-bg-prediction.pdf')