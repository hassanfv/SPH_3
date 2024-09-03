
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


f = h5py.File('grid_noneq_evolution_NeuralNet_T_6.2.hdf5', 'r')
TemperatureEvolution = f['TemperatureEvolution'][:]
AbundEvol = f['AbundanceEvolution'][:]

t_Arr_in_sec = f['TimeArray_seconds'][:]
t_Arr_in_yrs = t_Arr_in_sec / (3600. * 24. * 365.25)

N_nH = n_densities = f['TableBins/N_Densities'][()]

N_Z = n_metallicities = f['TableBins/N_Metallicities'][()]
N_T = n_temperatures = f['TableBins/N_Temperatures'][()]
densities = f['TableBins/Densities'][:]
metallicities = f['TableBins/Metallicities'][:]
temperatures = f['TableBins/Temperatures'][:]


nElm = 4

nXi = AbundEvol[0, 0, 0, nElm, :]

TEvol = TemperatureEvolution[0, 0, 0, :]
nHx = densities[0]
print('nH = ', 10**nHx)
plt.scatter(np.log10(TEvol), np.log10(nXi) , label = f'1', s = 5)



#===============================================================================================================


f = h5py.File('grid_noneq_evolution_NeuralNet_T_6.2x.hdf5', 'r')
TemperatureEvolution = f['TemperatureEvolution'][:]
AbundEvol = f['AbundanceEvolution'][:]

t_Arr_in_sec = f['TimeArray_seconds'][:]
t_Arr_in_yrs = t_Arr_in_sec / (3600. * 24. * 365.25)

N_nH = n_densities = f['TableBins/N_Densities'][()]

N_Z = n_metallicities = f['TableBins/N_Metallicities'][()]
N_T = n_temperatures = f['TableBins/N_Temperatures'][()]
densities = f['TableBins/Densities'][:]
metallicities = f['TableBins/Metallicities'][:]
temperatures = f['TableBins/Temperatures'][:]

nXi = AbundEvol[0, 0, 0, nElm, :]

TEvol = TemperatureEvolution[0, 0, 0, :]
nHx = densities[0]
print('nH = ', 10**nHx)
plt.scatter(np.log10(TEvol), np.log10(nXi) , label = f'2', s = 5)

#plt.yscale('log')

#plt.xlim(4, 10.5)
#plt.ylim(-11, -4)

plt.legend()

plt.savefig('Abund_vs_TEvol.png')


plt.show()








