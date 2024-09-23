
import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import h5py


#----- getMu
def getMu(Ab):

  s = 0.0
  p = 0.0
  for j in range(157):
    s += Ab[j] * AtomicMass[j]
    p += Ab[j] # Note that ne is also included in the sum!!

  mu = s / p
  
  return mu


df = pd.read_csv('data_species.csv')
print(df)
    
AtomicMass = df['A']

gamma = 5./3.
kB = 1.3807e-16
mH = 1.6726e-24
pc_to_cm = 3.086e18

iElm = 1 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

plt.figure(figsize = (8, 8))


#------ Single Chimes Output -----
f = h5py.File(f'SingleChimesRun.hdf5', 'r')
TempEvol = f['TemperatureEvolution'][:][0, 0, 0, :]
print('TempEvol.shape = ', TempEvol.shape)
AbundEvol = f['AbundanceEvolution'][:]     # (T, nH, Z, Elm, t)
print('AbundEvol.shape = ', AbundEvol.shape)
Abchim = AbundEvol[0, 0, 0, :, 0]
print(Abchim)
t_Arr_in_sec = f['TimeArray_seconds'][:]
t_Arr_in_yrs = t_Arr_in_sec / (3600. * 24. * 365.25)
#---------------------------------

AbundEvolX = AbundEvol[0, 0, 0, :, :]

AbEvol = AbundEvol[0, 0, 0, iElm, :]


#------ Finding the turning point -----
Tmp1 = TempEvol[:-1]
Tmp2 = TempEvol[1:]

nxx = np.where(((Tmp1 - Tmp2) < 0.0) & (Tmp2 < 1e4))[0] # Assuming turning point always occurs below 1e4 K. 

print(nxx)

if nxx.size > 0:
  ndx = nxx[0]
  T_crit = TempEvol[ndx]
  AbX = AbundEvolX[:, ndx]
  mu = getMu(AbX)
  print('ndx = ', ndx)
  print('T_crit = ', T_crit)
  print('mu_crit = ', mu)
else:
  print('nxx is empty !!!!!')


plt.scatter(np.log10(TempEvol), np.log10(1e-30+AbEvol), s = 5, color = 'b')

if nxx.size > 0:
  plt.axhline(y = np.log10(1e-30+AbEvol[ndx]), linestyle = '--', color = 'k')


plt.savefig('figx.png', dpi = 300, bbox_inches = 'tight') # Great! This confirms the idea !!!

plt.show()




