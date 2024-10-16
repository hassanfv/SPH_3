
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


df = pd.read_csv('xdata_species.xsv')
print(df)
    
AtomicMass = df['A']

gamma = 5./3.
kB = 1.3807e-16
mH = 1.6726e-24
pc_to_cm = 3.086e18


#------ Single Chimes Output -----
f = h5py.File(f'SingleChimesRun.hdf5', 'r')
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

ndx = 0

xRes = []

for ndx in range(0, AbEvol.shape[1], 1):

  AbX = AbEvol[:, ndx]
  mu = getMu(AbX)
  
  xRes.append([TempEvol[ndx], mu])

  print(f'T = {TempEvol[ndx]:.4E},   mu = {mu:.4f}')

xRes = np.array(xRes)

T = xRes[:, 0]
mu = xRes[:, 1]

n_steps = int((1.25 - 0.59) / (np.mean([0.04, 0.01])))
print('len n_steps = ', n_steps)
muSteps = np.linspace(0.04, 0.01, n_steps)


# Create a figure with 2 subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))

# First scatter plot
ax1.scatter(t_Arr_in_yrs, np.log10(TempEvol), color='k', s=20)
ax1.set_title("Plot 1: t_Arr_in_yrs vs log10(TempEvol)")
ax1.set_xlabel("t_Arr_in_yrs")
ax1.set_ylabel("log10(TempEvol)")
ax1.set_ylim(3.0, 10.5)

# Second scatter plot
ax2.scatter(np.log10(T), mu, color='k', s=10)
ax2.set_title("Plot 2: log10(T) vs mu")
ax2.set_xlabel("log10(T)")
ax2.set_ylabel("mu")
ax2.set_xlim(3.5, 10.5)
ax2.set_ylim(0.5, 1.3)

y0 = 0.59
i = 0
for tmp in muSteps:
  ax2.axhline(y = y0, linestyle = ':', color = 'b')
  y0 += tmp
  print(y0)
  i += 1

plt.savefig('figx.png', dpi = 300, bbox_inches = 'tight') # Great! This confirms the idea !!!

plt.show()




