
import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import time
import pickle
from numba import njit
import struct
import glob

'''
AbundanceEvolution
TableBins
TemperatureEvolution
TimeArray_seconds
'''


#f = h5py.File('./hdf5_files/grid_noneq_evolution_NeuralNet_rkpc_0.50_NH_19.4.hdf5', 'r')
f = h5py.File('grid_noneq_evolution_NeuralNetX_r_kpc_0.55_NH_21.25.hdf5', 'r')

# Print the attributes of HDF5 objects
for name, obj in f.items():
  print(name)
  for key, val in obj.attrs.items():
    print("    %s: %s" % (key, val))

print('\n\n')

 
TemperatureEvolution = f['TemperatureEvolution'][:]
print(f'TemperatureEvolution.shape = {TemperatureEvolution.shape}\n')
#print(f'TemperatureEvolution = {TemperatureEvolution}\n')

AbundanceEvolution = f['AbundanceEvolution'][:]
print('AbundanceEvolution.shape = ', AbundanceEvolution.shape)

with open('datax.pkl', 'wb') as filex:
  pickle.dump(TemperatureEvolution, filex)

t_Arr_in_sec = f['TimeArray_seconds'][:]
t_Arr_in_yrs = t_Arr_in_sec / (3600. * 24. * 365.25)

N_nH = n_densities = f['TableBins/N_Densities'][()]
print("N_Densities:", n_densities)

N_Z = n_metallicities = f['TableBins/N_Metallicities'][()]
print("N_Metallicities:", n_metallicities)

N_T = n_temperatures = f['TableBins/N_Temperatures'][()]
print("N_Temperatures:", n_temperatures)
print()
densities = f['TableBins/Densities'][:]
metallicities = f['TableBins/Metallicities'][:]
temperatures = f['TableBins/Temperatures'][:]

print('temperatures = ', temperatures)

inH = -1
print('lognH = ', densities[inH])
print('nH = ', 10**densities[inH])
iTemp = 0  #20 --> 1e6
print('T = ', 10**temperatures[iTemp])
iZ = 0
print('Z = ', metallicities[iZ])

#N_elem = 32 # 32 is from e upto Om i.e. from electron upto Neutral Oxygen

TEvol = TemperatureEvolution[iTemp, inH, iZ, :]


#--- This file can be used in prog1.py to compare the exact solution with the interpolated one! -------
tmp = {'t': t_Arr_in_yrs, 'T': np.log10(TEvol)}
with open('tmp.pkl', 'wb') as f:
  pickle.dump(tmp, f)
#------------------------------------------------------------------------------------------------------

plt.scatter(t_Arr_in_yrs, np.log10(TEvol), s = 5)

plt.ylim(3.75, 5.1)

plt.show()








