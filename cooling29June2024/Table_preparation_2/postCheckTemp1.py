
import numpy as np
#import pandas as pd
import pickle
import matplotlib.pyplot as plt
import h5py

pc_to_cm = 3.086e18

#------ Single Chimes Output -----
f = h5py.File(f'./SingleChimesRun.hdf5', 'r')
TempEvol = f['TemperatureEvolution'][:]
print(TempEvol.shape)
Tchim = np.log10(TempEvol[0, 0, 0, :])
print(Tchim)
#---------------------------------

plt.scatter(np.arange(101), Tchim, s = 5, color = 'k')
plt.ylim(3, 8)
plt.show()
s()


with open('TEvol.pkl', 'rb') as f:
  data = pickle.load(f)

#['TEvol', 'nH', 'Temp', 'rkpc', 'Lsh', 't']

#TEvol = np.zeros((N_rkpc, N_Lsh, N_nH, N_T, N_time))
#AbEvol = np.zeros((N_rkpc, N_Lsh, N_nH, N_T, N_Species, N_time))

TEvol = data['TEvol']
nH = data['nH']
Temp = data['Temp']
rkpc = data['rkpc']
Lsh = data['Lsh'] # ---- This is in log10 of parsec ! So 0 means 1.0 pc!
t = data['t']

print(TEvol.shape)

print()
print(f'nH = {nH}\n')
print(f'Temp = {Temp}\n')
print(f'rkpc = {rkpc}\n')
print(f'Lsh = {Lsh}\n')



# Helper function to find the bounding indices for interpolation
def find_bounds(arr, value):
    lower_idx = np.max(np.where(arr <= value))
    upper_idx = np.min(np.where(arr >= value))
    return lower_idx, upper_idx

# Linear interpolation between two points
def linear_interpolate(x0, x1, y0, y1, x):
    return y0 + (y1 - y0) * ((x - x0) / (x1 - x0))

# Multidimensional interpolation
def multidimensional_interpolate(rkpc, Lsh, nH, T, TEvol, rkpc_x, Lsh_x, nH_x, T_x):
    # Find the indices of the grid points that bound the target point in each dimension
    rkpc_lower, rkpc_upper = find_bounds(rkpc, rkpc_x)
    Lsh_lower, Lsh_upper = find_bounds(Lsh, Lsh_x)
    nH_lower, nH_upper = find_bounds(nH, nH_x)
    T_lower, T_upper = find_bounds(T, T_x)
    
    # Extract the TEvol arrays at the bounding grid points
    # Adjust the order of nH and Lsh as per the new specification
    TEvol0000 = TEvol[rkpc_lower, Lsh_lower, nH_lower, T_lower]
    TEvol0001 = TEvol[rkpc_lower, Lsh_lower, nH_lower, T_upper]
    TEvol0010 = TEvol[rkpc_lower, Lsh_lower, nH_upper, T_lower]
    TEvol0011 = TEvol[rkpc_lower, Lsh_lower, nH_upper, T_upper]
    TEvol0100 = TEvol[rkpc_lower, Lsh_upper, nH_lower, T_lower]
    TEvol0101 = TEvol[rkpc_lower, Lsh_upper, nH_lower, T_upper]
    TEvol0110 = TEvol[rkpc_lower, Lsh_upper, nH_upper, T_lower]
    TEvol0111 = TEvol[rkpc_lower, Lsh_upper, nH_upper, T_upper]
    TEvol1000 = TEvol[rkpc_upper, Lsh_lower, nH_lower, T_lower]
    TEvol1001 = TEvol[rkpc_upper, Lsh_lower, nH_lower, T_upper]
    TEvol1010 = TEvol[rkpc_upper, Lsh_lower, nH_upper, T_lower]
    TEvol1011 = TEvol[rkpc_upper, Lsh_lower, nH_upper, T_upper]
    TEvol1100 = TEvol[rkpc_upper, Lsh_upper, nH_lower, T_lower]
    TEvol1101 = TEvol[rkpc_upper, Lsh_upper, nH_lower, T_upper]
    TEvol1110 = TEvol[rkpc_upper, Lsh_upper, nH_upper, T_lower]
    TEvol1111 = TEvol[rkpc_upper, Lsh_upper, nH_upper, T_upper]
    
    # Perform interpolation in each dimension step by step
    TEvol00 = linear_interpolate(T[T_lower], T[T_upper], TEvol0000, TEvol0001, T_x)
    TEvol01 = linear_interpolate(T[T_lower], T[T_upper], TEvol0010, TEvol0011, T_x)
    TEvol10 = linear_interpolate(T[T_lower], T[T_upper], TEvol0100, TEvol0101, T_x)
    TEvol11 = linear_interpolate(T[T_lower], T[T_upper], TEvol0110, TEvol0111, T_x)
    TEvol0 = linear_interpolate(nH[nH_lower], nH[nH_upper], TEvol00, TEvol01, nH_x)
    TEvol1 = linear_interpolate(nH[nH_lower], nH[nH_upper], TEvol10, TEvol11, nH_x)
    TEvol_low = linear_interpolate(Lsh[Lsh_lower], Lsh[Lsh_upper], TEvol0, TEvol1, Lsh_x)
    
    TEvol00 = linear_interpolate(T[T_lower], T[T_upper], TEvol1000, TEvol1001, T_x)
    TEvol01 = linear_interpolate(T[T_lower], T[T_upper], TEvol1010, TEvol1011, T_x)
    TEvol10 = linear_interpolate(T[T_lower], T[T_upper], TEvol1100, TEvol1101, T_x)
    TEvol11 = linear_interpolate(T[T_lower], T[T_upper], TEvol1110, TEvol1111, T_x)
    TEvol0 = linear_interpolate(nH[nH_lower], nH[nH_upper], TEvol00, TEvol01, nH_x)
    TEvol1 = linear_interpolate(nH[nH_lower], nH[nH_upper], TEvol10, TEvol11, nH_x)
    TEvol_high = linear_interpolate(Lsh[Lsh_lower], Lsh[Lsh_upper], TEvol0, TEvol1, Lsh_x)
    
    TEvol_final = linear_interpolate(rkpc[rkpc_lower], rkpc[rkpc_upper], TEvol_low, TEvol_high, rkpc_x)
    
    return TEvol_final


rkpc_x = 0.25001
Lsh_x = 1.50001 # This is log10 of Lsh in pc!
nH_x = 3.0001
T_x = 5.60001
TEvol_x = multidimensional_interpolate(rkpc, Lsh, nH, Temp, TEvol, rkpc_x, Lsh_x, nH_x, T_x)

rkpc_x = 0.231
TEvol_x2= multidimensional_interpolate(rkpc, Lsh, nH, Temp, TEvol, rkpc_x, Lsh_x, nH_x, T_x)

rkpc_x = 0.281
TEvol_x3= multidimensional_interpolate(rkpc, Lsh, nH, Temp, TEvol, rkpc_x, Lsh_x, nH_x, T_x)


TT = TEvol[2, 6, 70, 26, :] # (N_rkpc, N_Lsh, N_nH, N_T, N_time)
TT2 = TEvol[3, 6, 70, 26, :] # (N_rkpc, N_Lsh, N_nH, N_T, N_time)

plt.scatter(t, TT, s = 5, color = 'k')
plt.scatter(t, TT2, s = 5, color = 'k')
plt.scatter(t, Tchim, s = 5, color = 'r')
plt.scatter(t, TEvol_x, s = 5, color = 'orange')
plt.scatter(t, TEvol_x2, s = 5, color = 'lime')
plt.scatter(t, TEvol_x3, s = 5, color = 'blue')

plt.savefig('fig2.png', bbox_inches = 'tight', dpi = 150)

plt.show()



