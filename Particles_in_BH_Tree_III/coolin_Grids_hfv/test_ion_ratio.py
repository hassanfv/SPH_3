
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

print()

inH = 0
print('nH = ', 10**densities[inH])
iTemp = 20 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

print()
print(f'nH0/nH = {(nH0[-1]/(nH0[-1]+nHp[-1])):.4f}, nHp/nH = {(nHp[-1]/(nH0[-1]+nHp[-1])):.4f}')
print(f'log(nH0/nH) = {np.log10(nH0[-1]/(nH0[-1]+nHp[-1])):.4f}, log(nHp/nH) = {np.log10(nHp[-1]/(nH0[-1]+nHp[-1])):.4f}')
print()
print(f'sort nelec = {np.sort(nelec)}')
print()


T = TEvol
nHeTot = nHe0 + nHep + nHepp

dictx = {'t_Arr_in_yrs': t_Arr_in_yrs, 'TEvol': TEvol, 'nHe0': nHe0, 'nHep': nHep, 'nHepp': nHepp, 'nH0': nH0, 'nHp': nHp,
         'nC0': nC0, 'nC1': nC1, 'nC2': nC2, 'nC3': nC3, 'nC4': nC4, 'nC5': nC5, 'nC6': nC6}

with open('chimesRes.pkl', 'wb') as f:
  pickle.dump(dictx, f)



#------ Result from "test_primordial_hdf5_v2.py" code -----
with open('chimesRes.pkl', 'rb') as f:
  df = pickle.load(f)
# dictx = {'t_Arr_in_yrs': t_Arr_in_yrs, 'TEvol': TEvol, 'nHe0': nHe0, 'nHep': nHep, 'nHepp': nHepp}
t_Arr_in_yrsx = df['t_Arr_in_yrs']
TEvolx = df['TEvol']
nH0x = df['nH0']
nHpx = df['nHp']
nHx = nH0x + nHpx
nHe0x = df['nHe0']
nHepx = df['nHep']
nHeppx = df['nHepp']
nHeTotx = nHe0x + nHepx + nHeppx

nC0x = df['nC0']
nC1x = df['nC1']
nC2x = df['nC2']
nC3x = df['nC3']
nC4x = df['nC4']
nC5x = df['nC5']
nC6x = df['nC6']
nCx = nC0x + nC1x + nC2x + nC3x + nC4x + nC5x + nC6x
#----------------------------------------------------------

nex = nHpx + (nHepx + 2.0 * nHeppx) + (nC1x + 2.0 * nC2x + 3.0 * nC3x + 4.0 * nC4x + 5.0 * nC5x + 6.0 * nC6x)

plt.plot(TEvolx, nC0x/nCx, color = 'r', )
plt.plot(TEvolx, nC1x/nCx, color = 'g', )
plt.plot(TEvolx, nC2x/nCx, color = 'b', )
plt.plot(TEvolx, nC3x/nCx, color = 'orange', )
plt.plot(TEvolx, nC4x/nCx, color = 'purple', )
plt.plot(TEvolx, nC5x/nCx, color = 'lime', )
plt.plot(TEvolx, nC6x/nCx, color = 'pink', )

plt.yscale('log')
plt.xscale('log')
plt.ylim(2e-3, 1.2)
plt.xlim(1e4, 1e6)

plt.title('Z = ' + str(metallicities[iZ]))

#plt.legend()

plt.tight_layout()

plt.savefig('C_ion_ratio.png')

plt.show()






