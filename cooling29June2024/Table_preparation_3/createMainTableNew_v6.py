
# In _New_v6 version, we use coarser grid for muG.
# In _New_v4 version, we are filling mainTable.
# In _New_v3 version, we include mu in the analysis.

import numpy as np
import matplotlib.pyplot as plt
import pickle
import glob
from scipy.interpolate import interp1d
import random
import pandas as pd
import os
import readchar

#===== closestNdx
def closestNdx(arr, val):
  return np.argmin(abs(arr - val))

#===== closestNdx
#def closestNdx(arr, val):
#  nx = np.argmin(abs(arr - val))
#  if arr[nx] > val: # We always get the left-hand side index!
#    nx = nx - 1
#  return nx


pc_to_cm = 3.086e18

#------------
#dfx = pd.read_csv('someInfoTable.csv')
#print(dfx.keys())
# ['nH', 'rkpc', 'Lsh', 't10', 't9', 't8', 't7', 't6', 't5', 't4', 't3', 'min_T', 'TR1', 'TR2']
#nHz = dfx['nH']
#rkpcz = dfx['rkpc']
#Lshz = dfx['Lsh']
#min_T = dfx['min_T']
#TR1 = dfx['TR1']
#TR2 = dfx['TR2']
#-------------


#----------- Preparing the grid -------
rkpcG = np.arange(0.01, 1.02, 0.1)# it is in kpc
LshG = np.arange(0.0, 2.51, 0.25) # it is in log10 of pc so 0.0 mean 1pc or 3.086e18 cm ---> We take Lsh in the range 1.0 pc up to ~300 pc.
nHG = np.arange(-4.0, 4.01, 0.1)
TG = np.arange(2.0, 10.31, 0.1)

N_rkpc = len(rkpcG)
N_Lsh = len(LshG)
N_nH = len(nHG)
N_T = len(TG)
#--------------------------------------

#-------- Creating mu grid ------------
#n_steps = int((1.25 - 0.59) / (np.mean([0.04, 0.01])))
#print('len n_steps = ', n_steps)
#muSteps = np.linspace(0.04, 0.01, n_steps)
#y0 = 0.59
#muG = np.zeros_like(muSteps)
#for i, tmp in enumerate(muSteps):
#  muG[i] = y0
#  y0 += tmp

#muG = np.arange(0.6, 1.26, 0.05)

muG = np.array([0.6, 0.7, 0.8, 0.9, 0.95, 1.0] + [float(_) for _ in np.arange(1.01, 1.24, 0.01)])

N_mu = len(muG)
#--------------------------------------

print(muG)


totTime = 50
tStep = 1
NtotTime = int(totTime/tStep)
print(f'NtotTime = {NtotTime}\n')

mainTable = np.zeros((N_nH, N_rkpc, N_Lsh, N_mu, N_T, NtotTime))
mainTableMu = np.zeros((N_nH, N_rkpc, N_Lsh, N_mu, N_T, NtotTime))

print(mainTable.shape)


print(f'N_rkpc * N_Lsh * N_nH * N_mu * NtotTime = {N_rkpc * N_Lsh * N_nH * N_mu * NtotTime}')

print()
print('muG = ', muG)
print('len(muG) = ', len(muG))



dirX = '/home/pc/Desktop/N_body_2024/SPH_3/cooling29June2024/Table_preparation_3/pklFromEC2/'
filez = glob.glob(dirX + '*.pkl')
print(len(filez))
print()

j = 5771#random.randint(0, 9800) #12    #11

nam = filez[j]

print(j, nam)

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

plt.scatter(t_Arr_in_yrs, TEvol, s = 1, color = 'k')
plt.show()




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

#------ For test --> You can remove it later !
T_p = 5.6
nx = closestNdx(TEvol, T_p)
print(f'nx = {nx}')
#------

#i = 24143

ndxTG_test = []

count = 0 # Only for test!
kk = 0 # Only for test!
csv_res = []

for i in range(1, len(TEvol)-1, 2):

  Ta = float(TEvol[i-1])
  Tb = float(TEvol[i])
  Tc = float(TEvol[i+1])
  
  mua = float(muEvol[i-1])
  mub = float(muEvol[i])
  muc = float(muEvol[i+1])

  ta = float(t_Arr_in_yrs[i-1])
  tb = float(t_Arr_in_yrs[i])
  tc = float(t_Arr_in_yrs[i+1])
  
  print('i = ', i)
  print('XX = ', Ta, Tb, Tc)
  print('YY = ', ta, tb, tc)

  print(f'BEFORE tZERO Subtraction = {ta, tb, tc}')

  tZeroPoint = ta

  ta -= tZeroPoint
  tb -= tZeroPoint
# tc will be corrected for tZeroPoint below in the else part of the if clause! 
  
  print(f'AFTER tZERO Subtraction = {ta, tb, tc}')

  #--- collecting those T from TG that are in the [Ta, Tc] interval -----
  Tres = []
  ndx_in_TG = []
  for ndxTG, Ti in enumerate(TG):
    #print('HH = ', Ta, Ti, Tc)
    if (Ta < Ti < Tc) | (Tc < Ti < Ta):
      print(f'\nINSIDE IF = {Ta, Ti, Tc}\n')
      Tres.append(float(Ti))
      ndx_in_TG.append(ndxTG)


  print('Tres = ', Tres)
  print()
  
  #if len(Tres) > 0:
  #  s()

  if Tres:
    print()
    print('--------------------')
    print(Ta, Tb, Tc)
    print(Tres)
    print()
    print('(ta, tb, tc, tZeroPoint) = ', ta, tb, tc, tZeroPoint)
    print()
    print('Tc ZZZZZ = ', Tc)

  #-- If Tres is not empty we will perfom the processing for each element in Tres -----
  if Tres:
    print('Tc BBBBBB = ', Tc)
    for ndxTG, T_p in zip(ndx_in_TG, Tres): # We do the following for every T_p in TG.
      
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
        
        muc = muEvol[i+1:i+1+k+1]
        muc = [float(x) for x in muc]
      else:
        print('HHEERREE', Tc, Tres, k)
        Tc = [ float(TEvol[i+1]) ]
        muc = [ float(muEvol[i+1]) ]
        tc = [t_Arr_in_yrs[i+1] - tZeroPoint]
      
      tarr = [ta, tb] + tc
      Tarr = [Ta, Tb] + Tc
      muarr = [mua, mub] + muc
    
      print('!!!!!!!!!!!!!!!!!!!!!')
      print(f'tarr = {tarr}\n')
      print(f'Tarr = {Tarr}\n')
      print(f'muarr = {muarr}\n')
      
      ndxTG_test.append(ndxTG) # !!!! Just for testing!

      # Now FineGriding:
      tFine = np.linspace(tarr[0], tarr[-1], 1000)
      T_interp = np.interp(tFine, tarr, Tarr)
      mu_interp = np.interp(tFine, tarr, muarr)
      
      nx = closestNdx(T_interp, T_p)
      
      nx_in_muG = closestNdx(muG, mu_interp[nx])
      mu_p = muG[nx_in_muG]
      print()
      print(f'muEvol[nx] = {muEvol[nx]},  mu_p = {mu_p}')
      print()
      
      Ttmp = T_interp[nx] # Only for test purposes!
      print('\nT_p, Ttmp = ', T_p, Ttmp)
      
      tarrX = tFine[nx:]
      TarrX = T_interp[nx:]
      muarrX= mu_interp[nx:]
      t100 = tarrX[0] + np.arange(0, NtotTime) # Note the usage of tStep as a multiplier !
      T100 = np.interp(t100, tarrX, TarrX)
      mu100 = np.interp(t100, tarrX, muarrX)

      mainTable[ndx_nH, ndx_rkpc, ndx_Lsh, nx_in_muG, ndxTG, :] = T100
      mainTableMu[ndx_nH, ndx_rkpc, ndx_Lsh, nx_in_muG, ndxTG, :] = mu100
      print()

      #print('Press a key\n')
      #kb = readchar.readkey()
 
      print(f'ndx_nH = {ndx_nH}, ndx_rkpc = {ndx_rkpc}, ndx_Lsh = {ndx_Lsh}, nx_in_muG = {nx_in_muG},  ndxTG = {ndxTG}')
      print(f'nHG[ndx_nH] = {nHG[ndx_nH]}, rkpcG[ndx_rkpc] = {rkpcG[ndx_rkpc]}, LshG[ndx_Lsh] = {LshG[ndx_Lsh]}, \
              muG[nx_in_muG] = {muG[nx_in_muG]}, TG[ndxTG]={TG[ndxTG]}')
      print()
      csv_res.append([T_interp[nx], TG[ndxTG], T100[0], T100[-1], mu_interp[nx], muG[nx_in_muG], mu100[0], mu100[-1]])
      
      #if True & (TG[ndxTG] < 6.8):
      if True:
        plt.scatter(t_Arr_in_yrs, TEvol, s = 1, color = 'grey')
        plt.scatter(np.array(tarr)+tZeroPoint, Tarr, color = 'blue')
        plt.scatter(np.array(tFine)+tZeroPoint, T_interp, s = 10, color = 'k')
        plt.scatter(np.array(t100)+tZeroPoint, T100, s = 1, color = 'lime')
        
        plt.axhline(y = T_p, linestyle = ':', color = 'red')
        
        plt.xlim(765337, 780000) # 766926
        plt.ylim(3.6, 6.05)
        
        ttmp = tFine[0] + tZeroPoint
        filename = f'test_{round(TG[ndxTG], 2)}_{mu_p:.6f}.jpg'
        plt.savefig(f'./pics/{filename}', bbox_inches='tight', dpi=300)
        
        #plt.show()
        plt.close()
        
        count +=1

#---- Saving T evolution
mainTable = mainTable.astype(np.float32)
with open('mainTable.pkl', 'wb') as f:
  pickle.dump(mainTable, f)

del mainTable
#---- Saving Mu evolution
mainTableMu = mainTableMu.astype(np.float32)
with open('mainTableMu.pkl', 'wb') as f:
  pickle.dump(mainTableMu, f)

print()
print('-----------')
print(ndx_nH, ndx_rkpc, ndx_Lsh, nx_in_muG)
print('-----------\n')

df = pd.DataFrame(csv_res)
df.columns = ['T_interp[nx]', 'TG[ndxTG]', 'T_start', 'T_End', 'mu_interp[nx]', 'muG[nx_in_muG]', 'mu_start', 'mu_End']
df.to_csv('dataTest.csv', index = False)
print(df.head(68))

#print()
#print('---------------------------------------------')
#for k in ndxTG_test:
#  print(k, TG[k])




