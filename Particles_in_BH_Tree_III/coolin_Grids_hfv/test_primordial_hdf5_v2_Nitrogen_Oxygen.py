
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


f = h5py.File('grid_noneq_evolution_primordialX.hdf5', 'r')

# Print the attributes of HDF5 objects
for name, obj in f.items():
  print(name)
  for key, val in obj.attrs.items():
    print("    %s: %s" % (key, val))
    
TemperatureEvolution = f['TemperatureEvolution'][:]
print(TemperatureEvolution.shape)

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
iTemp = 0
print('T = ', 10**temperatures[iTemp])
iZ = 0
print('Z = ', metallicities[iZ])


TEvol = TemperatureEvolution[iTemp, inH, iZ, :]

nelec = AbundanceEvolution[iTemp, inH, iZ, 0, :] * 10**densities[inH]

nH0 = AbundanceEvolution[iTemp, inH, iZ, 1, :] * 10**densities[inH]
nHp = AbundanceEvolution[iTemp, inH, iZ, 2, :] * 10**densities[inH]

nHe0 = AbundanceEvolution[iTemp, inH, iZ, 4, :] * 10**densities[inH]
nHep = AbundanceEvolution[iTemp, inH, iZ, 5, :] * 10**densities[inH]
nHepp= AbundanceEvolution[iTemp, inH, iZ, 6, :] * 10**densities[inH]

nC0= AbundanceEvolution[iTemp, inH, iZ, 7, :] * 10**densities[inH]
nC1= AbundanceEvolution[iTemp, inH, iZ, 8, :] * 10**densities[inH]
nC2= AbundanceEvolution[iTemp, inH, iZ, 9, :] * 10**densities[inH]
nC3= AbundanceEvolution[iTemp, inH, iZ, 10, :] * 10**densities[inH]
nC4= AbundanceEvolution[iTemp, inH, iZ, 11, :] * 10**densities[inH]
nC5= AbundanceEvolution[iTemp, inH, iZ, 12, :] * 10**densities[inH]
nC6= AbundanceEvolution[iTemp, inH, iZ, 13, :] * 10**densities[inH]
nCm= AbundanceEvolution[iTemp, inH, iZ, 14, :] * 10**densities[inH]

nN0= AbundanceEvolution[iTemp, inH, iZ, 15, :] * 10**densities[inH]
nN1= AbundanceEvolution[iTemp, inH, iZ, 16, :] * 10**densities[inH]
nN2= AbundanceEvolution[iTemp, inH, iZ, 17, :] * 10**densities[inH]
nN3= AbundanceEvolution[iTemp, inH, iZ, 18, :] * 10**densities[inH]
nN4= AbundanceEvolution[iTemp, inH, iZ, 19, :] * 10**densities[inH]
nN5= AbundanceEvolution[iTemp, inH, iZ, 20, :] * 10**densities[inH]
nN6= AbundanceEvolution[iTemp, inH, iZ, 21, :] * 10**densities[inH]
nN7= AbundanceEvolution[iTemp, inH, iZ, 22, :] * 10**densities[inH]

nO0= AbundanceEvolution[iTemp, inH, iZ, 23, :] * 10**densities[inH]
nO1= AbundanceEvolution[iTemp, inH, iZ, 24, :] * 10**densities[inH]
nO2= AbundanceEvolution[iTemp, inH, iZ, 25, :] * 10**densities[inH]
nO3= AbundanceEvolution[iTemp, inH, iZ, 26, :] * 10**densities[inH]
nO4= AbundanceEvolution[iTemp, inH, iZ, 27, :] * 10**densities[inH]
nO5= AbundanceEvolution[iTemp, inH, iZ, 28, :] * 10**densities[inH]
nO6= AbundanceEvolution[iTemp, inH, iZ, 29, :] * 10**densities[inH]
nO7= AbundanceEvolution[iTemp, inH, iZ, 30, :] * 10**densities[inH]
nO8= AbundanceEvolution[iTemp, inH, iZ, 31, :] * 10**densities[inH]
nOm= AbundanceEvolution[iTemp, inH, iZ, 32, :] * 10**densities[inH]

print()
print('nC0 + ... + nC6 = ', (nC0 + nC1 + nC2 + nC3 + nC4 + nC5 + nC6)[0])
print('expected nC_tot = ', 10**-3.61 * 10**densities[inH])
print()

print()
print(f'nH0/nH = {(nH0[-1]/(nH0[-1]+nHp[-1])):.4f}, nHp/nH = {(nHp[-1]/(nH0[-1]+nHp[-1])):.4f}')
print(f'log(nH0/nH) = {np.log10(nH0[-1]/(nH0[-1]+nHp[-1])):.4f}, log(nHp/nH) = {np.log10(nHp[-1]/(nH0[-1]+nHp[-1])):.4f}')
print()
print(f'sort nelec = {np.sort(nelec)}')
print()


T = TEvol
nHeTot = nHe0 + nHep + nHepp
plt.plot(T, nHe0/nHeTot, label = 'nHe0')
plt.plot(T, nHep/nHeTot, label = 'nHep')
plt.plot(T, nHepp/nHeTot,label = 'nHepp')
plt.yscale('log')
plt.xscale('log')
plt.ylim(2e-3, 1.2)
plt.xlim(1e4, 5e6)
plt.legend()
plt.show()



plt.figure(figsize = (12, 6))

plt.subplot(1, 2, 1)
plt.scatter(t_Arr_in_yrs, np.log10(TEvol), s = 5, color = 'k')
plt.xlim(0, 20000)
plt.ylim(3, 8)


plt.subplot(1, 2, 2)
plt.scatter(t_Arr_in_yrs, nHe0, s = 5, color = 'orange', label = 'nHe0')
plt.scatter(t_Arr_in_yrs, nHep, s = 5, color = 'lime', label = 'nHep')
plt.scatter(t_Arr_in_yrs, nHepp, s = 5, color = 'yellow', label = 'nHepp')
plt.yscale('log')
plt.legend()

dictx = {'t_Arr_in_yrs': t_Arr_in_yrs, 'TEvol': TEvol, 'nHe0': nHe0, 'nHep': nHep, 'nHepp': nHepp, 'nH0': nH0, 'nHp': nHp,
         'nC0': nC0, 'nC1': nC1, 'nC2': nC2, 'nC3': nC3, 'nC4': nC4, 'nC5': nC5, 'nC6': nC6, 'nCm': nCm,
         'nN0': nN0, 'nN1': nN1, 'nN2': nN2, 'nN3': nN3, 'nN4': nN4, 'nN5': nN5, 'nN6': nN6, 'nN7': nN7,
         'nO0': nO0, 'nO1': nO1, 'nO2': nO2, 'nO3': nO3, 'nO4': nO4, 'nO5': nO5, 'nO6': nO6, 'nO7': nO7, 'nO8': nO8, 'nOm': nOm}

dicLOG= {'t_Arr_in_yrs': t_Arr_in_yrs, 'TEvol': np.log10(TEvol), 'nHe0': np.log10(nHe0+1e-30), 'nHep': np.log10(nHep+1e-30),
         'nHepp': np.log10(nHepp+1e-30), 'nH0': np.log10(nH0+1e-30), 'nHp': np.log10(nHp+1e-30), 'nC0': np.log10(nC0+1e-30),
         'nC1': np.log10(nC1+1e-30), 'nC2': np.log10(nC2+1e-30), 'nC3': np.log10(nC3+1e-30), 'nC4': np.log10(nC4+1e-30),
         'nC5': np.log10(nC5+1e-30), 'nC6': np.log10(nC6+1e-30)}

with open('chimesRes_C_N_O.pkl', 'wb') as f:
  pickle.dump(dictx, f)

#with open('chimesResLOG.pkl', 'wb') as f:
#  pickle.dump(dicLOG, f)

plt.savefig('primordial.png')

plt.show()






