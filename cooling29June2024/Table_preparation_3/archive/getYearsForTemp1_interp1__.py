
import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import h5py
import time


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


TTA = time.time()

#------ Single Chimes Output -----
#f = h5py.File(f'SingleChimesRun.hdf5', 'r')
f = h5py.File(f'./nH_2.4_rkpc_0.81_Lsh_20.739.hdf5', 'r')
TempEvol = f['TemperatureEvolution'][:][0, 0, 0, :]
print('TempEvol.shape = ', TempEvol.shape)
AbundEvol = f['AbundanceEvolution'][:]     # (T, nH, Z, Elm, t)
print('AbundEvol.shape = ', AbundEvol.shape)
Abchim = AbundEvol[0, 0, 0, :, 0]
#print(Abchim)
t_Arr_in_sec = f['TimeArray_seconds'][:]
t_Arr_in_yrs = t_Arr_in_sec / (3600. * 24. * 365.25)
#---------------------------------

AbEvol = AbundEvol[0, 0, 0, :, :]
print(AbEvol.shape)


Res = np.zeros(9)

checker1 = checker2 = checker3 = checker4 = checker5 = checker6 = checker7 = checker8 = 0

#=========================
for i, T in enumerate(TempEvol):
  if (np.log10(T) < 10.0) and (checker1 == 0):
    Res[0] = t_Arr_in_yrs[i]
    checker1 = 1
  
  if (np.log10(T) < 9.0) and (checker2 == 0):
    Res[1] = t_Arr_in_yrs[i]
    checker2 = 1
    
  if (np.log10(T) < 8.0) and (checker3 == 0):
    Res[2] = t_Arr_in_yrs[i]
    checker3 = 1
  
  if (np.log10(T) < 7.0) and (checker4 == 0):
    Res[3] = t_Arr_in_yrs[i]
    checker4 = 1
    
  if (np.log10(T) < 6.0) and (checker5 == 0):
    Res[4] = t_Arr_in_yrs[i]
    checker5 = 1
  
  if (np.log10(T) < 5.0) and (checker6 == 0):
    Res[5] = t_Arr_in_yrs[i]
    checker6 = 1
    
  if (np.log10(T) < 4.0) and (checker7 == 0):
    Res[6] = t_Arr_in_yrs[i]
    checker7 = 1
  
  if (np.log10(T) < 3.0) and (checker8 == 0):
    Res[7] = t_Arr_in_yrs[i]
    checker8 = 1
#=========================


ndx = 0
xRes = []

for ndx in range(0, AbEvol.shape[1], 1):

  AbX = AbEvol[:, ndx]
  mu = getMu(AbX)
  
  xRes.append([TempEvol[ndx], mu])

  #print(f'T = {TempEvol[ndx]:.4E},   mu = {mu:.4f}')

xRes = np.array(xRes)

T = xRes[:, 0]
mu = xRes[:, 1]

n_steps = int((1.25 - 0.59) / (np.mean([0.04, 0.01])))
print('len n_steps = ', n_steps)
muSteps = np.linspace(0.04, 0.01, n_steps)

print(f'Years in yrs for different T = {Res} ')

print()
print(f'Elapsed time = {time.time() - TTA:.2f} seconds.')
print()

# Create a figure with 2 subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))

# First scatter plot
ax1.scatter(t_Arr_in_yrs, np.log10(TempEvol), color='k', s=20)
ax1.set_title("Plot 1: t_Arr_in_yrs vs log10(TempEvol)")
ax1.set_xlabel("t_Arr_in_yrs")
ax1.set_ylabel("log10(TempEvol)")
ax1.set_ylim(2.0, 10.5)

for timex in Res:
  if timex > 0.0:
    ax1.axvline(x = timex, linestyle = ':', color = 'b')


# Second scatter plot
ax2.scatter(np.log10(T), mu, color='k', s=10)
ax2.set_title("Plot 2: log10(T) vs mu")
ax2.set_xlabel("log10(T)")
ax2.set_ylabel("mu")
ax2.set_xlim(2.0, 10.5)
ax2.set_ylim(0.5, 1.3)

y0 = 0.59
i = 0
for tmp in muSteps:
  ax2.axhline(y = y0, linestyle = ':', color = 'b')
  y0 += tmp
  #print(y0)
  i += 1

plt.savefig('figx.png', dpi = 300, bbox_inches = 'tight') # Great! This confirms the idea !!!

plt.show()




