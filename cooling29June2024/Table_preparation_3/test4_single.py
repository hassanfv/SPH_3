
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle
import glob


#----------- Preparing the grid -------
rkpcG = np.arange(0.01, 1.02, 0.1)# it is in kpc
LshG = np.arange(0.0, 2.51, 0.25) # it is in log10 of pc so 0.0 mean 1pc or 3.086e18 cm ---> We take Lsh in the range 1.0 pc up to ~300 pc.
nHG = np.arange(-4.0, 4.01, 0.1)
TG = np.arange(2.0, 10.31, 0.1)
#-------- Creating mu grid ------------
muG = np.arange(0.6, 1.21, 0.1)
#--------------------------------------

print(muG)



N_nH = len(nHG)
N_rkpc = len(rkpcG)
N_Lsh = len(LshG)
N_T = len(TG)
N_mu = len(muG)

#------ TEST FILE ----------------------
dirX = '/home/pc/Desktop/N_body_2024/SPH_3/cooling29June2024/Table_preparation_3/pklFromEC2/'
filez = glob.glob(dirX + '*.pkl')
print(len(filez))
print()
j = 11
nam = filez[j]
with open(nam, 'rb') as f:
  data = pickle.load(f)
nH_i = nH_p = float(data['nH_p'])
rkpc_i = rkpc_p = float(data['rkpc_p'])
Lsh_i = Lsh_p = float(data['Lsh_p'])
print(f'nH_p = {nH_p},  rkpc_p = {rkpc_p},  Lsh_p = {Lsh_p}\n')

TEvol = np.log10(data['TempEvol'])
t_Arr_in_yrs = data['t_in_sec'] / 3600. / 24. / 365.25
#---------------------------------------

with open('mainTable.pkl', 'rb') as f:
  data = pickle.load(f)

with open('mainTableMu.pkl', 'rb') as f:
  dataMu = pickle.load(f)
# mainTable = np.zeros((N_nH, N_rkpc, N_Lsh, N_mu, N_T, 101))

#--- T_i = 7.6 is our starting point
T_i = 3.7373 # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
mu_i = 1.04867 # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Lsh_i = 0.0
# Step 1: Find the index fo, nHi, rkpc_i, Lsh_i, mu_i, and T_i in nHG, rkpcG, LshG, muG, and TG !
# NOTE: We also need the evolved mu, i.e. the mu after the dt timestep evolution !

nx_nH_i = round((nH_i - nHG[0]) / 0.1)
nx_rkpc_i = round((rkpc_i - rkpcG[0]) / 0.1)
nx_Lsh_i = round((Lsh_i - LshG[0]) / 0.25)
nx_mu_i = np.argmin(abs(muG - mu_i))
nx_T_i = round((T_i - TG[0]) / 0.1)


k = 0
TEvolved = data[nx_nH_i, nx_rkpc_i, nx_Lsh_i, nx_mu_i+k, nx_T_i]
while TEvolved[0] == 0:
  
  print('mu in LOOP = ', muG[nx_mu_i+k])
  print('TEvolved in LOOP = ', TEvolved)
  TEvolved = data[nx_nH_i, nx_rkpc_i, nx_Lsh_i, nx_mu_i+k, nx_T_i]#[:20]   #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  MuEvolved = dataMu[nx_nH_i, nx_rkpc_i, nx_Lsh_i, nx_mu_i+k, nx_T_i]#[:20]#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  k += 1

delta_T = TEvolved[-1] - TEvolved[0]

delta_Mu= MuEvolved[-1] - MuEvolved[0]

print(f'T Before = {T_i}')
print(f'mu Before = {mu_i}')

print('TEvolved = ', TEvolved)






