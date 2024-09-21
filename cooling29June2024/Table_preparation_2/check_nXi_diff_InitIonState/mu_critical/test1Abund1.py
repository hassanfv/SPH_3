
import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import h5py


iElm = 10 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

AbEvol = AbundEvol[0, 0, 0, iElm, :]

plt.scatter(np.log10(TempEvol), np.log10(1e-30+AbEvol), s = 5, color = 'b')


plt.savefig('figx.png', dpi = 300, bbox_inches = 'tight') # Great! This confirms the idea !!!

plt.show()




