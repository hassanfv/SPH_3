
# NO NEED FOR MPI !!!!!!!

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import time
import pickle
import struct
import glob


filz = np.sort(glob.glob('./coolHeatGrids_pkl/*.pkl'))

print(filz)

kpc = np.array([np.float32(fil[-10:-7]) for fil in filz])

print(kpc)

j = 0 # We read the first file so that we can get N_ of different varaiables.

with open(filz[j], 'rb') as f:
  data = pickle.load(f)

#print(data.keys())
#['densities', 'metallicities', 'temperatures', 'timeArr_in_sec', 'uEvolution', 'muArr']

densities = data['densities']
metallicities = data['metallicities']
temperatures = data['temperatures']
timeArr_in_sec = data['timeArr_in_sec']
uEvolution = data['uEvolution']
muArr = data['muArr']
metalz = data['metalz']

print('uEvolution.shape = ', uEvolution.shape)
print('muArr.shape = ', muArr.shape)
print('metalz.shape = ', metalz.shape)

N_kpc = len(kpc)
N_T = len(temperatures)
N_nH = len(densities)
N_Z = len(metallicities)
N_M = len(metalz[0, 0, 0, :, 0])
N_time = len(timeArr_in_sec)

print('(N_kpc, N_T, N_nH, N_Z, N_M, N_time = )', N_kpc, N_T, N_nH, N_Z, N_M, N_time)

uEvolutionX = np.zeros((N_kpc, N_T, N_nH, N_Z, N_time))
muArrX = np.zeros((N_kpc, N_T, N_nH, N_Z, N_time))
metalzX= np.zeros((N_kpc, N_T, N_nH, N_Z, N_M, N_time))
print('uEvolutionX.shape = ', uEvolutionX.shape)
print('muArrX.shape = ', muArrX.shape)
print('metalzX.shape = ', metalzX.shape)


for j in range(N_kpc):

  with open(filz[j], 'rb') as f:
    data = pickle.load(f)

  uEvolution = data['uEvolution']
  muArr = data['muArr']
  metalz= data['metalz']
  
  uEvolutionX[j, :, :, :, :] = uEvolution
  muArrX[j, :, :, :, :] = muArr
  metalzX[j, :, :, :, :, :] = metalz

print()
for i in range(N_kpc):
  print(f'test_{i} = {uEvolutionX[i, 8, 20, 1, 2]}')
  print(f'test_{i} = {metalzX[i, 8, 20, 1, 1, 2]}')
  print(f'test_{i} = {metalzX[i, 8, 20, 1, 2, 2]}')
  print()
print()


#dictx = {'densities': densities, 'metallicities': metallicities, 'temperatures': temperatures, 'timeArr_in_sec': timeArr_in_sec,
#         'uEvolution': uEvolutionX, 'muArr': muArrX, 'kpc': kpc}
#with open('coolHeatGridX.pkl', 'wb') as f:
#  pickle.dump(dictx, f)



with open('coolHeatGridZ.bin', 'wb') as f:
  # Save N_Kpc, N_nH, N_Z, N_T, N_time for dimensions
  f.write(struct.pack('i', N_kpc))
  f.write(struct.pack('i', N_nH))
  f.write(struct.pack('i', N_Z))
  f.write(struct.pack('i', N_T))
  f.write(struct.pack('i', N_M))
  f.write(struct.pack('i', N_time)) # Note that len(muArr) = N_time!
  
  # Save kpc
  for val in kpc:
    f.write(struct.pack('f', val))  # use 'f' format for float
  
  # Save densities
  for val in densities:
    f.write(struct.pack('f', val))  # use 'f' format for float
      
  # Save metallicities
  for val in metallicities:
    f.write(struct.pack('f', val))  # use 'f' format for float
      
  # Save temperatures
  for val in temperatures:
    f.write(struct.pack('f', val))  # use 'f' format for float
      
  # Save timeArr
  for val in timeArr_in_sec:
    f.write(struct.pack('f', val))  # use 'f' format for float
      
  # Save res
  for t in range(N_kpc):
    for i in range(N_T):
      for j in range(N_nH):
        for k in range(N_Z):
          for l in range(N_time):
            f.write(struct.pack('f', uEvolutionX[t][i][j][k][l]))  # use 'f' format for float

  # Save muA
  for t in range(N_kpc):
    for i in range(N_T):
      for j in range(N_nH):
        for k in range(N_Z):
          for l in range(N_time):
            f.write(struct.pack('f', muArrX[t][i][j][k][l]))  # use 'f' format for float

  # Save metalz
  for t in range(N_kpc):
    for i in range(N_T):
      for j in range(N_nH):
        for k in range(N_Z):
          for ii in range(N_M):
            for l in range(N_time):
              f.write(struct.pack('f', metalzX[t][i][j][k][ii][l]))  # use 'f' format for float




