
import numpy as np
#import pandas as pd
import pickle
import matplotlib.pyplot as plt
import h5py
import glob
import time

pc_to_cm = 3.086e18

'''
#------ Single Chimes Output -----
f = h5py.File(f'./SingleChimesRun.hdf5', 'r')
TempEvol = f['TemperatureEvolution'][:]
print(TempEvol.shape)
Tchim = np.log10(TempEvol[0, 0, 0, :])
print(Tchim)
#---------------------------------
'''


pklFiles = glob.glob('./on_single_EC2/Out_pkl/*.pkl')

with open(pklFiles[0], 'rb') as f:
  data = pickle.load(f)
#['TempEvol', 'AbundEvol', 'nH', 'Temp', 'rkpc', 'Lsh', 'Species_id', 'Species_name', 't']

nHG = data['nH']
TempG = data['Temp']
rkpcG = data['rkpc']
LshG = data['Lsh'] # ---- This is in log10 of parsec ! So 0 means 1.0 pc!
tG = data['t']
Species_id = data['Species_id']
Species_name = data['Species_name']

TempEvol = data['TempEvol']
print('TempEvol.shape = ', TempEvol.shape)
print()
AbundEvol = data['AbundEvol']
print('AbundEvol.shape = ', AbundEvol.shape)

print()
print(f'nH = {nHG}\n')
print(f'Temp = {TempG}\n')
print(f'rkpc = {rkpcG}\n')
print(f'Lsh = {LshG}\n')
print(f'Species_id = {Species_id}\n')
print(f'Species_name = {Species_name}\n')


TA = time.time()

N_Species = len(Species_id)
N_rkpc = len(rkpcG)
N_Lsh = len(LshG)
N_T = len(TempG)
N_nH = len(nHG)
N_time = len(tG)

TEvol = np.zeros((N_rkpc, N_Lsh, N_nH, N_T, N_time))
AbEvol = np.zeros((N_rkpc, N_Lsh, N_nH, N_T, N_Species, N_time))

print(f'N_rkpc = {N_rkpc},   N_Lsh = {N_Lsh}')
print('N_rkpc * N_Lsh = ', N_rkpc * N_Lsh)

print()


for i in range(N_rkpc):
  for j in range(N_Lsh):
  
  
    file_pkl = f'./on_single_EC2/Out_pkl/rkpc_{rkpcG[i]:.2f}_Lsh_{np.log10(10**LshG[j] * pc_to_cm):.3f}.pkl'
    with open(file_pkl, 'rb') as f:
      data = pickle.load(f)
    
    TempEvol = data['TempEvol']
    AbundEvol = data['AbundEvol']
    
    print()
    print('TempEvol.shape = ', TempEvol.shape)
    print()
    
    TB = time.time()
    for k in range(N_nH):
      for l in range(N_T):
        TEvol[i, j, k, l, :] = np.log10(TempEvol[l, k, 0, :]) # TEvol[rkpc, Lsh, nH, T, time]
        for b in range(N_Species): # Note that AbundEvol contains only 12 elements instead of full 157 ! Selection was done in progSep2024_v8.py code!
          tmp = np.log10(1e-30+AbundEvol[l, k, b, :]) # Note that the Z dimension was already removed from AbundEvol in progSep2024_v8.py code!
          nt = np.where(tmp < -30.0)[0]
          if len(nt) > 0:
            tmp[nt] = -30.0
          AbEvol[i, j, k, l, b, :] = tmp


TEvolx = {'TEvol': TEvol, 'nH': nHG, 'Temp': TempG, 'rkpc': rkpcG, 'Lsh': LshG, 't': np.arange(N_time)}
with open('TEvol.pkl', 'wb') as f:
  pickle.dump(TEvolx, f)

AbEvolx = {'AbEvol': AbEvol, 'nH': nHG, 'Temp': TempG, 'rkpc': rkpcG, 'Lsh': LshG, 'Species_id': Species_id, 'Species_name': Species_name, 't': np.arange(N_time)}
with open('AbEvol.pkl', 'wb') as f:
  pickle.dump(AbEvolx, f)


print(f'Elapsed timeX = {(time.time() - TA):.1f} sec')


