
import numpy as np
#import pandas as pd
import pickle
import matplotlib.pyplot as plt
import h5py


iElm = 0 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

plt.figure(figsize = (8, 8))

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

AbEvol = AbundEvol[0, 0, 0, iElm, :]

plt.scatter(np.log10(TempEvol), np.log10(1e-30+AbEvol), s = 5, color = 'k')


#------ Single Chimes Output -----
f = h5py.File(f'./SingleChimesRunX.hdf5', 'r')
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

plt.scatter(np.log10(TempEvol), np.log10(1e-30+AbEvol), s = 5, color = 'b')


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

plt.scatter(np.log10(TempEvol), np.log10(1e-30+AbEvol), s = 5, color = 'lime')

'''
#------ Single Chimes Output -----
f = h5py.File(f'./SingleChimesRun_T_5.5_Ion_1.hdf5', 'r')
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

plt.scatter(np.log10(TempEvol), np.log10(1e-30+AbEvol), s = 5, color = 'gold')


#------ Single Chimes Output -----
f = h5py.File(f'./SingleChimesRun_T_5.5_Ion_5.hdf5', 'r')
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

plt.scatter(np.log10(TempEvol), np.log10(1e-30+AbEvol), s = 5, color = 'pink')


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

AbEvol = AbundEvol[0, 0, 0, iElm, :]

plt.scatter(np.log10(TempEvol), np.log10(1e-30+AbEvol), s = 5, color = 'red')

print('t_Arr_in_yrs = ', t_Arr_in_yrs)


#------ Single Chimes Output -----
f = h5py.File(f'./SingleChimesRun_T_6.0_Ion_5_nH_3.5.hdf5', 'r')
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

plt.scatter(np.log10(TempEvol), np.log10(1e-30+AbEvol), s = 15, color = 'purple')
'''

#plt.ylim(-8, 1)

plt.savefig('fig.png', dpi = 300, bbox_inches = 'tight') # Great! This confirms the idea !!!

plt.show()




