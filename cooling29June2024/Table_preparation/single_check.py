
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


f = h5py.File(f'./grid_noneq_evolution_NeuralNet.hdf5', 'r')

 
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

inH = 0
print('nH = ', 10**densities[inH])

iTemp = 0  #20 --> 1e6
print('T = ', 10**temperatures[iTemp])
iZ = 0
print('Z = ', metallicities[iZ])

TEvol = TemperatureEvolution[iTemp, inH, iZ, :]

plt.scatter(t_Arr_in_yrs, np.log10(TEvol), s = 5, color = 'k')


plt.ylim(3, 11)

#plt.xlim(0, 3000)

#plt.yscale('log')

plt.savefig('r_0.51.png')

plt.show()








