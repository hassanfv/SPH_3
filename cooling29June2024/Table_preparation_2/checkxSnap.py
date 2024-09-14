
import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import time
import pickle
import struct
import glob

'''
TempEvol.shape = (1517, 101)
AbundEvol.shape =  (1517, 157, 101)
'''


f = h5py.File(f'./grid_noneq_evolution_Snap.hdf5', 'r')

 
TempEvol = f['TemperatureEvolution'][:]
print(f'TempEvol.shape = {TempEvol.shape}\n')

AbundEvol = f['AbundanceEvolution']
print('AbundEvol.shape = ', AbundEvol.shape)


TEvol = TempEvol[171, :]

print()
print('TEvol = ', TEvol)

plt.scatter(np.arange(101), np.log10(TEvol), s = 5, color = 'k')

plt.ylim(3, 8)

plt.show()





