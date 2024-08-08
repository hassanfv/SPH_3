
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

rkpci = 0.8
NHj = 18

f = h5py.File(f'./hdf5_files/grid_noneq_evolution_NeuralNet_rkpc_{rkpci}_NH_{NHj}.hdf5', 'r')

# Print the attributes of HDF5 objects
for name, obj in f.items():
  print(name)
  for key, val in obj.attrs.items():
    print("    %s: %s" % (key, val))

print('\n\n')

 
TemperatureEvolution = f['TemperatureEvolution'][:]
print(f'TemperatureEvolution.shape = {TemperatureEvolution.shape}\n')
#print(f'TemperatureEvolution = {TemperatureEvolution}\n')

AbundanceEvolution = f['AbundanceEvolution']
print('AbundanceEvolution.shape = ', AbundanceEvolution.shape)

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

inH = 38
inH2= 39
print('nH = ', 10**densities[inH])
print('nH2 = ', 10**densities[inH2])
nHRatio = 10**densities[inH2] / 10**densities[inH]
print()
print(f'nHRatio = {nHRatio}')
print()

iTemp = 0  #20 --> 1e6
iTemp2 = 0  #20 --> 1e6
print('T = ', 10**temperatures[iTemp])
iZ = 0
print('Z = ', metallicities[iZ])

TEvol = TemperatureEvolution[iTemp, inH, iZ, :]
TEvol2= TemperatureEvolution[iTemp, inH2, iZ, :]


t_Arr_in_yrs2 = nHRatio * t_Arr_in_yrs

plt.scatter(t_Arr_in_yrs, TEvol, s = 5, color = 'k')
plt.scatter(t_Arr_in_yrs2, TEvol2, s = 5, color = 'orange')



plt.ylim(1e3, 5e5)

#plt.xlim(0, 3000)

plt.yscale('log')

plt.show()








