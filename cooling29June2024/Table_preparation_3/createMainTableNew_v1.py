
# In _v3 version, we check whether a T_p is within accepted range. someInfoTable.csv contains values like min_T, TR1 and TR2 than is used for this purpose.

import numpy as np
import matplotlib.pyplot as plt
import pickle
import glob
from scipy.interpolate import interp1d
import random
import pandas as pd

#===== closestNdx
def closestNdx(arr, val):
  return np.argmin(abs(arr - val))


#===== closestNdx
def boundNdx(arr, val):
  nx = np.argmin(abs(arr - val))
  if val >= arr[nx]:
    nxLeft = nx
    nxRight= nx + 1
  else:
    nxLeft = nx - 1
    nxRight = nx
  return nxLeft, nxRight


pc_to_cm = 3.086e18

#------------
dfx = pd.read_csv('someInfoTable.csv')
print(dfx.keys())
# ['nH', 'rkpc', 'Lsh', 't10', 't9', 't8', 't7', 't6', 't5', 't4', 't3', 'min_T', 'TR1', 'TR2']
nHz = dfx['nH']
rkpcz = dfx['rkpc']
Lshz = dfx['Lsh']
min_T = dfx['min_T']
TR1 = dfx['TR1']
TR2 = dfx['TR2']
#-------------


#----------- Preparing the grid -------
rkpcG = np.arange(0.01, 1.02, 0.1)# it is in kpc
LshG = np.arange(0.0, 2.51, 0.25) # it is in log10 of pc so 0.0 mean 1pc or 3.086e18 cm ---> We take Lsh in the range 1.0 pc up to ~300 pc.
nHG = np.arange(-4.0, 4.01, 0.1)
TG = np.arange(2.0, 10.31, 0.1)
#--------------------------------------

#-------- Creating mu grid ------------
n_steps = int((1.25 - 0.59) / (np.mean([0.04, 0.01])))
print('len n_steps = ', n_steps)
muSteps = np.linspace(0.04, 0.01, n_steps)
y0 = 0.59
muG = np.zeros_like(muSteps)
for i, tmp in enumerate(muSteps):
  muG[i] = y0
  y0 += tmp
#--------------------------------------

print()
print('muG = ', muG)
print('len(muG) = ', len(muG))



dirX = '/home/pc/Desktop/N_body_2024/SPH_3/cooling29June2024/Table_preparation_3/pklFromEC2/'
filez = glob.glob(dirX + '*.pkl')
print(len(filez))
print()

j = 12#random.randint(0, 9800) #12    #11

nam = filez[j]

print(nam)

with open(nam, 'rb') as f:
  data = pickle.load(f)
# 'TempEvol', 'AbundEvol', 'nH', 'rkpc', 'Lsh', 'Species_id', 'Species_name', 't_in_sec', 'nH_p', 'rkpc_p', 'Lsh_p', 'mu'
keyz = list(data.keys())

TEvol = np.log10(data['TempEvol'])
AbEvol = data['AbundEvol']
Species_id = data['Species_id']
Species_name = data['Species_name']
nH_p = float(data['nH_p'])
rkpc_p = float(data['rkpc_p'])
Lsh_p = float(data['Lsh_p'])
t_Arr_in_yrs = data['t_in_sec'] / 3600. / 24. / 365.25

print()
print(f'nH_p, rkpc_p, Lsh_p = {nH_p}, {rkpc_p}, {Lsh_p}')

print()
print(keyz)
print()

if 'mu' in keyz:
  muEvol = data['mu']
else:
  muEvol = 0.6 + np.zeros_like(TEvol)

print('muEvol = ', muEvol)
print()

#----- Finding the indices for nH, rkpc, Lsh, and mu to be used in mainTable -------
ndx_nH = closestNdx(nHG, nH_p)
ndx_rkpc = closestNdx(rkpcG, rkpc_p)
ndx_Lsh = closestNdx(LshG, Lsh_p)
#-----------------------------------------------------------------------------------

print('TEvol = ', TEvol)

#------ For test
T_p = 5.6
nx = closestNdx(TEvol, T_p)
print(f'nx = {nx}')
#------

i = 24143

#for i in range(1, len(TEvol)-1):

Ta = float(TEvol[i-1])
Tb = float(TEvol[i])
Tc = float(TEvol[i+1])

ta = float(t_Arr_in_yrs[i-1])
tb = float(t_Arr_in_yrs[i])
tc = float(t_Arr_in_yrs[i+1])

tZeroPoint = ta

ta -= tZeroPoint
tb -= tZeroPoint

#--- collecting those T from TG that are in the [Ta, Tc] interval -----
Tres = []
for Ti in TG:
  if (Ta <= Ti <= Tc) | (Tc <= Ti <= Ta):
    Tres.append(float(Ti))


if Tres:
  print()
  print('--------------------')
  print(Ta, Tb, Tc)
  print(Tres)
  print()
  print('(ta, tb, tc, tZeroPoint) = ', ta, tb, tc, tZeroPoint)
  print()

#-- If Tres is not empty we will perfom the processing for each element in Tres -----
if Tres:
  for T_p in Tres: # We do the following for every T_p in TG.
    
    dt = 0.0
    k = 0
    while dt < 500:
      k += 1
      dt = t_Arr_in_yrs[i+1+k] - t_Arr_in_yrs[i+1]
    
    if k > 1:
      tc = t_Arr_in_yrs[i+1:i+1+k+1]
      tc = [float(x) for x in tc]
      tc = [x - tZeroPoint for x in tc]
      
      Tc = TEvol[i+1:i+1+k+1]
      Tc = [float(x) for x in Tc]
    else:
      Tc = [ float(Tc) ] # tZeroPoint is subtracted!
    
    tarr = [ta, tb] + tc
    Tarr = [Ta, Tb] + Tc
  
    print('!!!!!!!!!!!!!!!!!!!!!')
    print(f'tarr = {tarr}')
    print(f'Tarr = {Tarr}')

    # Now FineGriding:
    tFine = np.linspace(tarr[0], tarr[-1], 1000)
    T_interp = np.interp(tFine, tarr, Tarr)
    
    nx = closestNdx(T_interp, T_p)
    
    Ttmp = T_interp[nx] # Only for test purposes!
    print('\nT_p, Ttmp = ', T_p, Ttmp)
    
    tarrX = tFine[nx:]
    TarrX = T_interp[nx:]
    t100 = tarrX[0] + np.arange(0, 101)
    T100 = np.interp(t100, tarrX, TarrX)
    
    plt.scatter(tarr, Tarr, color = 'blue')
    plt.scatter(tFine, T_interp, s = 10, color = 'k')
    plt.scatter(t100, T100, s = 1, color = 'lime')
    
    plt.axhline(y = T_p, linestyle = ':', color = 'red')
    plt.show()







