
import numpy as np
#import pandas as pd
import pickle
import matplotlib.pyplot as plt
import h5py


#iElm = 3 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#------ Single Chimes Output -----
f = h5py.File(f'./SingleChimesRun_T_6.0_Ion_5.hdf5', 'r')
TempEvol = f['TemperatureEvolution'][:]
print('TempEvol.shape = ', TempEvol.shape)
AbundEvol = f['AbundanceEvolution'][:]     # (T, nH, Z, Elm, t)
print('AbundEvol.shape = ', AbundEvol.shape)
Abchim = AbundEvol[0, 0, 0, :, 0]
print(Abchim)
t_Arr_in_sec = f['TimeArray_seconds'][:]
t_Arr_in_yrs = t_Arr_in_sec / (3600. * 24. * 365.25)
#---------------------------------

plt.scatter(t_Arr_in_yrs, np.log10(TempEvol), s = 5, color = 'k')

plt.show()




