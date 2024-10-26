
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

j = random.randint(0, 9800) #12    #11

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

T_p = 5.6

nx = closestNdx(TEvol, T_p)

dt = t_Arr_in_yrs[nx+1] - t_Arr_in_yrs[nx]
print(f'\nt1, t2, t3 before adding k (nx + k) = {0, t_Arr_in_yrs[nx]-t_Arr_in_yrs[nx-1], t_Arr_in_yrs[nx+1]-t_Arr_in_yrs[nx-1]}')
print(f'T1, T2, T3 before adding k (nx + k) = {TEvol[nx-1], TEvol[nx], TEvol[nx+1]}\n')
print(f'dt before: {dt}\n')
k = 1
while dt < 500: # We want to make sure we have enough data points after nx (at least 100 years... but we assume 500 to be safe!)
  dt = t_Arr_in_yrs[nx+k] - t_Arr_in_yrs[nx]
  k += 1

print(f'dt after: {dt}\n')
print(f'k = {k}\n')

t1 = float(t_Arr_in_yrs[nx-1])
tZeroPoint = t1+0.0
t2 = float(t_Arr_in_yrs[nx])
if k > 1:
  t3 = t_Arr_in_yrs[nx+1:nx+k+1]
  t3 = [float(x) for x in t3]
  t3 = [x - tZeroPoint for x in t3]
else:
  t3 = [ float(t_Arr_in_yrs[nx+1]) - tZeroPoint ] # tZeroPoint is subtracted!

t1 -= tZeroPoint
t2 -= tZeroPoint
tarr = [t1, t2] + t3

T1 = float(TEvol[nx-1])
T2 = float(TEvol[nx])

#----
if k > 1: # we do this because TEvol[nx+1:nx+1] is EMPTY but TEvol[nx+1] is not Empty!
  T3 = TEvol[nx+1:nx+k+1]
  T3 = [float(x) for x in T3]
else:
  T3 = [ float(TEvol[nx+1]) ]
#----

Tarr = [T1, T2] + T3

print(f'\ntarr Original = {tarr}\n')
print(f'\nTarr Original = {Tarr}\n')
print(f'Sizes before ploting: len(tarr), len(Tarr) = {len(tarr), len(Tarr)}')

plt.scatter(tarr, Tarr, color = 'blue')

delta_T = np.abs(T_p - T2)

print(f'\ndelta_T = {delta_T}\n')
print('\n---------------')
print(f'First closest T to {T_p} is {T2}\n')
print(f'\nXXX: T_p = {T_p},  T[nx] = {TEvol[nx]}\n')


if delta_T < 0.01:
  print('\n!!!!!!!!!!!!! delta_T < 0.01 !!!!!!!!!!!!!!')
  t100 = tarr[1] + np.arange(0, 101)
  T100 = np.interp(t100, tarr[1:], Tarr[1:])

if delta_T >= 0.01: # We fine grid the Tarr so that we have a value that is close enough to T_p ---> Good approach !
  print('\n!!!!!!!!!!!!! delta_T > 0.01 !!!!!!!!!!!!!!')
  tFine = np.linspace(tarr[0], tarr[-1], 1000)
  T_interp = np.interp(tFine, tarr, Tarr)
  
  nx = closestNdx(T_interp, T_p)
  
  T2 = T_interp[nx]

  delta_T = np.abs(T_p - T2)
  print(f'\ndelta_T = {delta_T}\n\n')
  
  print(f'Second closest T to {T_p} is {T2}\n')
  print()
  print(f'tFine[0] = {tFine[0]}, tFine[-1] = {tFine[-1]}')
  print()
  
  tarrX = tFine[nx:]
  TarrX = T_interp[nx:]
  t100 = tarrX[0] + np.arange(0, 101)
  T100 = np.interp(t100, tarrX, TarrX)

  print(f'\nYYY: T_p = {T_p},  T[nx] = {T_interp[nx]}\n')
  
  plt.scatter(tFine, T_interp, s = 10, color = 'k')

plt.scatter(t100, T100, s = 1, color = 'lime')
plt.axhline(y = T_p, linestyle = ':', color = 'red')
plt.show()






