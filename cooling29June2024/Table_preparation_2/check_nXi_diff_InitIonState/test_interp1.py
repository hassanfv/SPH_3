
import numpy as np
#import pandas as pd
import pickle
import matplotlib.pyplot as plt
import h5py
from scipy.interpolate import interp1d


iElm = 1 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#------ Single Chimes Output -----
f = h5py.File(f'./SingleChimesRun.hdf5', 'r')
TempEvol = f['TemperatureEvolution'][:]
print('TempEvol.shape = ', TempEvol.shape)
AbundEvol = f['AbundanceEvolution'][:]     # (T, nH, Z, Elm, t)
print('AbundEvol.shape = ', AbundEvol.shape)
Abchim = AbundEvol[0, 0, 0, :, 0]
print(Abchim)
t_Arr_in_sec = f['TimeArray_seconds'][:]
t_Arr_in_yrs = t_Arr_in_sec / (3600. * 24. * 365.25)
#---------------------------------

AbEvol = AbundEvol[0, 0, 0, iElm, :]

T = np.log10(TempEvol[0, 0, 0, :])
A = np.log10(1e-30+AbEvol)

print(T.shape, A.shape)

Tgrid = np.linspace(3.0, 10.2, 1000)
interp_func = interp1d(T, A, bounds_error=False, fill_value=(A[0], A[-1]))
A2 = interp_func(Tgrid)

plt.scatter(T, A, s = 5, color = 'k')
#plt.scatter(Tgrid, A2, s = 5, color = 'blue')


plt.savefig('figX.png', dpi = 300, bbox_inches = 'tight') # Great! This confirms the idea !!!

plt.show()




