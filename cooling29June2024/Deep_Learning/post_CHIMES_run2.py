
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
iTemp = 20
print('T = ', 10**temperatures[iTemp])
iZ = 0
print('Z = ', metallicities[iZ])


#TEvol = TemperatureEvolution[iTemp, inH, iZ, :]

N_elem = 32 # 32 is from e upto Om i.e. from electron upto Neutral Oxygen

Xres = []

T1 = time.time()

for inH in range(N_nH):

  for iTemp in range(N_T):

    for j in range(len(t_Arr_in_yrs)):

      TEvol = TemperatureEvolution[iTemp, inH, iZ, :]
      TEvol = [np.log10(_) for _ in TEvol]

      res_tmp = [t_Arr_in_yrs[j], densities[inH], TEvol[j]]
      
      for i in range(N_elem):
      
        Abundx = AbundanceEvolution[iTemp, inH, iZ, i, j] * 10**densities[inH]
        
        if Abundx == 0.0:
          Abundx = -30.0
        else:
          Abundx = np.log10(Abundx)
        
        res_tmp += [Abundx]

      Xres.append(res_tmp)

Xres = np.array(Xres)

print(Xres)
print()
print(Xres.shape)

with open('initialData.pkl', 'wb') as f:
  pickle.dump(Xres, f)

print('Elapsed time = ', time.time() - T1)





