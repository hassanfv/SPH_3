
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


#----- plot_h
def plot_h(inH):
  TEvol = TemperatureEvolution[0, inH, 0, :]
  nHx = 10**densities[inH]
  print('nH = ', nHx)
  plt.plot(t_Arr_in_yrs, np.log10(TEvol), label = f'nH = {nHx:.1f} cm^-3', linewidth = 2)
  plt.ylim(3.9, 5.1)


plot_h(75)
plot_h(76)
plot_h(77)
plot_h(78)
plot_h(79)
plot_h(80)

plt.xlabel('t (years)')
plt.ylabel('logT')

plt.legend()

plt.show()



