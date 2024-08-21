
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

def find_closest(T, T0):

    differences = np.abs(T - T0)
    index = np.argmin(differences)
    closest_value = T[index]

    return index, closest_value


f = h5py.File('grid_noneq_evolution_NeuralNetX_r_kpc_0.5.hdf5', 'r')
TemperatureEvolution = f['TemperatureEvolution'][:]
AbundEvol = f['AbundanceEvolution'][:]

print('Abund.shape = ', AbundEvol.shape)

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

print()
#print('TEvol = ', TEvol)

T0 = 250000.0
ndx_t, Tx = find_closest(TEvol, T0)
print(Tx)

nElm = 5
nXi = AbundEvol[0, 0, 0, nElm, ndx_t]
print('nXi = ', nXi)


plt.scatter(np.log10(TEvol), AbundEvol[0, 0, 0, nElm, :], label = f'rkpc = 0.5', s = 5)


#plt.scatter(t_Arr_in_yrs, np.log10(TEvol), label = f'rkpc = 0.5', s = 5)
#plt.legend()
#plt.savefig('diff_rkpc_same_T.png')

plt.show()








