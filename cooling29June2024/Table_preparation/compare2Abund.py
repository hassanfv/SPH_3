
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
f = h5py.File('grid_noneq_evolution_NeuralNet.hdf5', 'r')

# Print the attributes of HDF5 objects
for name, obj in f.items():
  print(name)
  for key, val in obj.attrs.items():
    print("    %s: %s" % (key, val))

print('\n\n')

 
TemperatureEvolution = f['TemperatureEvolution'][:]
print(f'TemperatureEvolution.shape = {TemperatureEvolution.shape}\n')
#print(f'TemperatureEvolution = {TemperatureEvolution}\n')

#with open('datax.pkl', 'wb') as f:
#  pickle.dump(TemperatureEvolution, f)
#s()

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


#----- plot_h
def plot_h1(inH, iElm):
  AbundEvol = AbundanceEvolution[0, inH, 0, iElm, :] # 11 is O+4 or OV.
  nHx = densities[inH]
  print('nH = ', 10**nHx)
  #plt.plot(t_Arr_in_yrs*10, np.log10(AbundEvol), label = f'lognH = {nHx:.2f} cm^-3', linewidth = 2, color = 'b')
  plt.scatter(t_Arr_in_yrs, np.log10(AbundEvol), label = f'lognH = {nHx:.2f} cm^-3', color = 'b', s = 5)


#----- plot_h
def plot_h2(inH, iElm):
  AbundEvol = AbundanceEvolution[0, inH, 0, iElm, :] # 11 is O+4 or OV.
  nHx = densities[inH]
  print('nH = ', 10**nHx)
  #plt.plot(t_Arr_in_yrs, np.log10(AbundEvol), label = f'lognH = {nHx:.2f} cm^-3', linewidth = 2, color = 'orange')
  plt.scatter(t_Arr_in_yrs, np.log10(AbundEvol), label = f'lognH = {nHx:.2f} cm^-3', color = 'orange', s = 5)
  

iElm = 13

plot_h1(0, iElm)
plot_h2(0, iElm)

#plt.xlim(0, 1e6)

plt.legend()

plt.show()








