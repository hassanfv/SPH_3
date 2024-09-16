
import numpy as np
#import pandas as pd
import pickle
import matplotlib.pyplot as plt
import h5py


#iElm = 3 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

plt.scatter(t_Arr_in_yrs, np.log10(TempEvol), s = 5, color = 'k')


plt.ylim(3.7, 6.2)
plt.xlim(0.0, 8000)

plt.xlabel('Time [yrs]')
plt.ylabel('log(T)')

plt.savefig('T_vs_time.png', dpi = 300, bbox_inches = 'tight') # Great! This confirms the idea !!!



plt.show()




