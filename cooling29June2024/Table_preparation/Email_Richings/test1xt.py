
import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import time
import pickle
from numba import njit
import struct
import glob

f = h5py.File('grid_noneq_evolution_NeuralNet.hdf5', 'r')
 
TemperatureEvolution = f['TemperatureEvolution'][:]

t_Arr_in_sec = f['TimeArray_seconds'][:]
t_Arr_in_yrs = t_Arr_in_sec / (3600. * 24. * 365.25)

densities = f['TableBins/Densities'][:]
metallicities = f['TableBins/Metallicities'][:]
temperatures = f['TableBins/Temperatures'][:]


#----- plot_h
def plot_h(inH):
  TEvol = TemperatureEvolution[0, inH, 0, :]
  nHx = densities[inH]
  print('nH = ', nHx)
  #plt.plot(t_Arr_in_yrs, np.log10(TEvol), label = f'lognH = {nHx:.2f} cm^-3', linewidth = 2)
  plt.scatter(t_Arr_in_yrs, np.log10(TEvol), label = f'lognH = {nHx:.2f} cm^-3', s = 1)


plot_h(5)
plot_h(-1)

plt.xlabel('t (years)')
plt.ylabel('logT')

plt.legend()

plt.savefig('plotx.png', dpi = 200, bbox_inches = 'tight')

plt.show()



