import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import h5py
import os

pc_to_cm = 3.086e18

#----------- Preparing the grid -------
rkpcG = np.arange(0.01, 1.02, 0.1)# it is in kpc
LshG = np.arange(0.0, 2.51, 0.25) # it is in log10 of pc so 0.0 mean 1pc or 3.086e18 cm ---> We take Lsh in the range 1.0 pc up to ~300 pc.
nHG = np.arange(-4.0, 4.01, 0.1)
#--------------------------------------
print(len(rkpcG), len(LshG), len(nHG), len(rkpcG)*len(LshG)*len(nHG))
N_nH = len(nHG)
N_rkpc = len(rkpcG)
N_Lsh = len(LshG)
inList = [] #np.zeros(N_nH * N_rkpc * N_Lsh)

for i in range(N_nH):
  for j in range(N_rkpc):
    for k in range(N_Lsh):
      inList.append([float(nHG[i]), float(rkpcG[j]), float(LshG[k])])

N = len(inList)

#ResX = np.zeros((N, 2))

for j in range(N):

  nH_t, rkpc_t, Lsh_t = inList[j]

  nam = f'./nH_{nH_t:.1f}_rkpc_{rkpc_t:.2f}_Lsh_{np.log10(10**Lsh_t * pc_to_cm):.3f}.pkl'
  nam = 'nH_0.8_rkpc_0.61_Lsh_20.239.pkl'
  
  print(nam)
  
  with open(nam, 'rb') as f:
    data = pickle.load(f)
  # ['TempEvol', 'AbundEvol', 'nH', 'rkpc', 'Lsh', 'Species_id', 'Species_name', 't_in_sec']
    
  TempEvol = data['TempEvol']
  AbundEvol = data['AbundEvol']
  nH = data['nH']
  rkpc = data['rkpc']
  Lsh = data['Lsh']
  Species_id = data['Species_id']
  Species_name = data['Species_name']
  t_in_sec = data['t_in_sec']
  
  print(f'TempEvol.shape = {TempEvol.shape}')
  print(f'AbundEvol.shape = {AbundEvol.shape}')
  
  t_in_yrs = t_in_sec / 3600.0 / 24.0 / 365.25
  
  plt.scatter(t_in_yrs, np.log10(TempEvol), color='k', s=10)
  plt.ylim(2.0, 10.5)
  
  plt.show()
  
  s()




