
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

 
TempEvol = f['TemperatureEvolution'][:]
print(f'TempEvol.shape = {TempEvol.shape}\n')

AbundEvol = f['AbundanceEvolution']
print('AbundEvol.shape = ', AbundEvol.shape)

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
print()
print('densities = ', densities)

'''
TempEvol.shape = (37, 41, 1, 101)
AbundEvol.shape =  (37, 41, 1, 157, 101)
N_Densities: 41
N_Metallicities: 1
N_Temperatures: 37
'''

ndx_T = 4
ndx_nH = 39

TEvol = TempEvol[ndx_T, ndx_nH, 0, :]

print()
print('TEvol = ', TEvol)

plt.scatter(t_Arr_in_yrs, np.log10(TEvol), s = 5, color = 'k', label = f'nH = {(10**densities[ndx_nH]):.2f}')
plt.legend()

plt.show()





