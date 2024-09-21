
import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import h5py

#----- Temp_to_u
def Temp_to_u(T, Ab):

  s = 0.0
  p = 0.0
  for j in range(157):
    s += Ab[j] * AtomicMass[j]
    p += Ab[j] # Note that ne is also included in the sum!!

  mu = s / p

  utmp = kB / mH / (gamma - 1.) / mu * T
  
  return utmp, mu


df = pd.read_csv('data_species.csv')
print(df)
    
AtomicMass = df['A']

gamma = 5./3.
kB = 1.3807e-16
mH = 1.6726e-24




iElm = 1 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

plt.figure(figsize = (8, 8))


#------ Single Chimes Output -----
f = h5py.File(f'./SingleChimesRun_initIonState_6.hdf5', 'r')
TempEvol = f['TemperatureEvolution'][:][0, 0, 0, :]
print('TempEvol.shape = ', TempEvol.shape)
AbundEvol = f['AbundanceEvolution'][:]     # (T, nH, Z, Elm, t)
print('AbundEvol.shape = ', AbundEvol.shape)
Abchim = AbundEvol[0, 0, 0, :, 0]
print(Abchim)
t_Arr_in_sec = f['TimeArray_seconds'][:]
t_Arr_in_yrs = t_Arr_in_sec / (3600. * 24. * 365.25)
#---------------------------------

AbEvol = AbundEvol[0, 0, 0, iElm, :]

print(TempEvol)
print()


#------ Finding the turning point -----
for j in range(len(TempEvol)-1):
  dT = TempEvol[j] - TempEvol[j+1]
  #print(i, TempEvol[i], TempEvol[i+1], dT)
  if dT < 0.0:
    ndx = j
    print('ndx = ', j)
    break

uX = np.zeros_like(TempEvol)
muX = np.zeros_like(TempEvol)

for i in range(0, len(TempEvol)):
  
  AbX = AbundEvol[0, 0, 0, :, i]
  TX = TempEvol[i]
  u, mu = Temp_to_u(TX, AbX)
  
  #print(f'{TX}, {mu:.3f}, {u:.3E}')
  uX[i] = u
  muX[i] = mu

print('logu = ', np.log10(uX))


#plt.scatter(np.log10(uX), np.log10(1e-30+AbEvol), s = 5, color = 'b')

plt.scatter(np.log10(TempEvol), muX, s = 5, color = 'b')

plt.axhline(y = muX[ndx], linestyle = '--', color = 'k')



plt.savefig('fig_uThermal.png', dpi = 300, bbox_inches = 'tight') # Great! This confirms the idea !!!

plt.show()




