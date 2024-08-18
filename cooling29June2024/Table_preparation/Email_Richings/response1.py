
import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import time
import pickle
from numba import njit
import struct
import glob

f = h5py.File('grid_noneq_evolution_NeuralNet_rkpc_0.60.hdf5', 'r')
 
TemperatureEvolution = f['TemperatureEvolution'][:]

t_Arr_in_sec = f['TimeArray_seconds'][:]
t_Arr_in_yrs = t_Arr_in_sec / (3600. * 24. * 365.25)

densities = f['TableBins/Densities'][:]
metallicities = f['TableBins/Metallicities'][:]
temperatures = f['TableBins/Temperatures'][:]

AbundanceEvolution = f['AbundanceEvolution']
print('AbundanceEvolution.shape = ', AbundanceEvolution.shape) #-----> (1, 81, 1, 157, 5001)--> (ndx_T, ndx_nH, ndx_Z, ndx_elm, ndx_time)


ndx_nH = 80
ndx_elm = 1 # 1 ---> HI fraction, i.e. HI/Htot


def plot_HIfrac(ndx_nH, ndx_elm):
  nHx = densities[ndx_nH]
  xAbund = AbundanceEvolution[0, ndx_nH, 0, ndx_elm, :]
  plt.plot(t_Arr_in_yrs, np.log10(xAbund), label = f'lognH = {nHx:.2f} cm^-3', linewidth = 2)
  #plt.ylim(0.0, 1.1)


#----- plot_h
def plot_h(inH):
  TEvol = TemperatureEvolution[0, inH, 0, :]
  nHx = densities[inH]
  print('nH = ', nHx)
  plt.plot(t_Arr_in_yrs, np.log10(TEvol), label = f'lognH = {nHx:.2f} cm^-3', linewidth = 0.5)
  plt.ylim(3.9, 5.1)



plot_HIfrac(75, 1)
plot_HIfrac(76, 1)
plot_HIfrac(77, 1)
plot_HIfrac(78, 1)
plot_HIfrac(79, 1)
plot_HIfrac(80, 1)

plt.xlabel('t (years)')
plt.ylabel('HI fraction')

plt.legend()

plt.savefig('plotHI.png', dpi = 200, bbox_inches = 'tight')

plt.show()



