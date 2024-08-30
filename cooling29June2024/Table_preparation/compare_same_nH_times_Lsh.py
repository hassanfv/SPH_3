
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


f = h5py.File('grid_noneq_evolution_NeuralNet_nH_2_Lsh_17.hdf5', 'r')
TemperatureEvolution = f['TemperatureEvolution'][:]
AbundanceEvolution = f['AbundanceEvolution'][:]

t_Arr_in_sec = f['TimeArray_seconds'][:]
t_Arr_in_yrs = t_Arr_in_sec / (3600. * 24. * 365.25)

N_nH = n_densities = f['TableBins/N_Densities'][()]

N_Z = n_metallicities = f['TableBins/N_Metallicities'][()]
N_T = n_temperatures = f['TableBins/N_Temperatures'][()]
densities = f['TableBins/Densities'][:]
metallicities = f['TableBins/Metallicities'][:]
temperatures = f['TableBins/Temperatures'][:]

inH = 0

TEvol = TemperatureEvolution[0, inH, 0, :]
nHx = densities[inH]
print('nH = ', 10**nHx)
plt.scatter(t_Arr_in_yrs, np.log10(TEvol), label = f'nH_2_Lsh_17', s = 5)



#===============================================================================================================
f = h5py.File('grid_noneq_evolution_NeuralNet_nH_3_Lsh_16.hdf5', 'r')
TemperatureEvolution = f['TemperatureEvolution'][:]
AbundanceEvolution = f['AbundanceEvolution'][:]

t_Arr_in_sec = f['TimeArray_seconds'][:]
t_Arr_in_yrs = t_Arr_in_sec / (3600. * 24. * 365.25)

N_nH = n_densities = f['TableBins/N_Densities'][()]

N_Z = n_metallicities = f['TableBins/N_Metallicities'][()]
N_T = n_temperatures = f['TableBins/N_Temperatures'][()]
densities = f['TableBins/Densities'][:]
metallicities = f['TableBins/Metallicities'][:]
temperatures = f['TableBins/Temperatures'][:]

inH = 0

TEvol = TemperatureEvolution[0, inH, 0, :]
nHx = densities[inH]
print('nH = ', 10**nHx)
plt.scatter(t_Arr_in_yrs, np.log10(TEvol), label = f'nH_3_Lsh_16', s = 5)



plt.legend()

plt.savefig('Same_nH_times_Lsh.png')

plt.show()








