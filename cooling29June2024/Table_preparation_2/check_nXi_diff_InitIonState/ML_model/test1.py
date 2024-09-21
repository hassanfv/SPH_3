
import numpy as np
#import pandas as pd
import pickle
import matplotlib.pyplot as plt
import h5py
from scipy.interpolate import interp1d


iElm = 1 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#------ Single Chimes Output -----
f = h5py.File(f'./SingleChimesRun.hdf5', 'r')
TempEvol = f['TemperatureEvolution'][:][0, 0, 0, :]
print('TempEvol.shape = ', TempEvol.shape)
AbundEvol = f['AbundanceEvolution'][:]     # (T, nH, Z, Elm, t)
print('AbundEvol.shape = ', AbundEvol.shape)


SpName = ['HI', 'OI', 'CI', 'CII', 'CIV', 'SiII', 'SiIV', 'NV', 'FeII', 'MgI', 'MgII', 'SII']
#                HI OI CI CII CIV SiII SiIV  NV  FeII  MgI  MgII  SII
SelectSpecies = [1, 23, 7, 8, 10, 58,  60,   19, 111,  44,  45,   73]

AbundEvol = f['AbundanceEvolution'][:][0, 0, 0, SelectSpecies, :]     # (T, nH, Z, Elm, t)
print('AbundEvol.shape = ', AbundEvol.shape)


N = len(TempEvol)
X = np.zeros((N, 13))
Y = np.zeros((N, 13))

for i in range(N-1):
  X[i, :12] = np.log10(1e-30+AbundEvol[:, i])
  X[i, 12] = TempEvol[i]
  
  Y[i, :12] = np.log10(1e-30+AbundEvol[:, i+1])
  Y[i, 12] = TempEvol[i+1]

print(X)
print()
print(Y)


from sklearn.model_selection import train_test_split
from sklearn.neural_network import MLPRegressor
from sklearn.metrics import mean_squared_error

X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.2, random_state=42)
model = MLPRegressor(hidden_layer_sizes=(64, 32), activation='relu', solver='adam', max_iter=1000, random_state=42)

# Train the model
model.fit(X_train, y_train)
y_pred = model.predict(X_test)

mse = mean_squared_error(y_test, y_pred)
print(f"Mean Squared Error: {mse}")

print(y_pred.shape) # (2001, 13)

plt.scatter(y_test[:, 1], y_pred[:, 1])

plt.show()







