# Pol Benítez Colominas, July 2024
# Universitat Politècnica de Catalunya

# Script to train basic machine learning models and predict band gaps for CAP solid solutions

import numpy as np
import pandas as pd

from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.linear_model import LinearRegression
from sklearn.kernel_ridge import KernelRidge
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score

import joblib

# Import the data
df_solid_solutions = pd.read_csv('bg-solid-solutions.csv')

S_conc = np.array(df_solid_solutions['S_concentration'])
Br_conc = np.array(df_solid_solutions['Br_concentration'])
bg = np.array(df_solid_solutions['bandgap'])

# Define X, and y tensors
X = []
for x in range(len(S_conc)):
    row = [S_conc[x], Br_conc[x]]
    X.append(row)
X = np.array(X)

y = bg

# Split X, y in train and test set
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.30, random_state=42)

##########  Random Forest Regressor
# Train the model
model = RandomForestRegressor()

model.fit(X_train, y_train)

# Make predictions
y_pred = model.predict(X_test)

# Evaluate the model
mse = mean_squared_error(y_test, y_pred)
mae = mean_absolute_error(y_test, y_pred)
r2 = r2_score(y_test, y_pred)

print(f'Random Forest Regressor\nMSE: {mse}\nMAE: {mae}\nr^2: {r2}\n')
##########################################################

##########  Gradient Boosting Regressor
# Train the model
model = GradientBoostingRegressor()

model.fit(X_train, y_train)

# Make predictions
y_pred = model.predict(X_test)

# Evaluate the model
mse = mean_squared_error(y_test, y_pred)
mae = mean_absolute_error(y_test, y_pred)
r2 = r2_score(y_test, y_pred)

print(f'Gradient Boosting Regressor\nMSE: {mse}\nMAE: {mae}\nr^2: {r2}\n')
##########################################################

##########  Linear Regression
# Train the model
model = LinearRegression()

model.fit(X_train, y_train)

# Make predictions
y_pred = model.predict(X_test)

# Evaluate the model
mse = mean_squared_error(y_test, y_pred)
mae = mean_absolute_error(y_test, y_pred)
r2 = r2_score(y_test, y_pred)

print(f'Linear Regression\nMSE: {mse}\nMAE: {mae}\nr^2: {r2}\n')
##########################################################

##########  Kernel Ridge
# Train the model
model = KernelRidge(kernel='rbf')

model.fit(X_train, y_train)

# Make predictions
y_pred = model.predict(X_test)

# Evaluate the model
mse = mean_squared_error(y_test, y_pred)
mae = mean_absolute_error(y_test, y_pred)
r2 = r2_score(y_test, y_pred)

print(f'Kernel Ridge\nMSE: {mse}\nMAE: {mae}\nr^2: {r2}\n')
##########################################################

##########  Gaussian Process Regressor
# Train the model
kernel = C(1.0, (1e-3, 1e3)) * RBF(1.0, (1e-3, 1e3))
model = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=10, random_state=42)

model.fit(X_train, y_train)

# Make predictions
y_pred = model.predict(X_test)

# Evaluate the model
mse = mean_squared_error(y_test, y_pred)
mae = mean_absolute_error(y_test, y_pred)
r2 = r2_score(y_test, y_pred)

print(f'Gaussian Process Regressor\nMSE: {mse}\nMAE: {mae}\nr^2: {r2}\n')
##########################################################

# Save the model
joblib.dump(model, 'trained_model.joblib')