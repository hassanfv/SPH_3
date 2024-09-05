
import numpy as np
#import pandas as pd
import pickle
import matplotlib.pyplot as plt
import h5py


with open('AbEvol.pkl', 'rb') as f:
  data = pickle.load(f)

#['AbEvol', 'nH', 'Temp', 'rkpc', 'Lsh', 'Species_id', 'Species_name', 't']

#TEvol = np.zeros((N_rkpc, N_Lsh, N_nH, N_T, N_time))
#AbEvol = np.zeros((N_rkpc, N_Lsh, N_nH, N_T, N_Species, N_time))

AbEvol = 10**data['AbEvol'] # (6, 6, 41, 37, 12, 101)
nH = 10**data['nH']
Temp = 10**data['Temp']
rkpc = data['rkpc']
Lsh = data['Lsh'] # ---- This is in log10 of parsec ! So 0 means 1.0 pc!
t = data['t']

Species_id = data['Species_id']
Species_name = data['Species_name']
print(AbEvol.shape)

iElm = 0 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#------ Single Chimes Output -----
f = h5py.File(f'./SingleChimesRun.hdf5', 'r')
#TempEvol = f['TemperatureEvolution'][:]
AbundEvol = f['AbundanceEvolution'][:]     # (T, nH, Z, Elm, t)
print(AbundEvol.shape)
Abchim = np.log10(AbundEvol[0, 0, 0, Species_id[iElm], :])
print(Abchim)
#---------------------------------

print()
print(f'nH = {nH}\n')
print(f'Temp = {Temp}\n')
print(f'rkpc = {rkpc}\n')
print(f'Lsh = {Lsh}\n')
print(f'Species_id = {Species_id}\n')
print(f'Species_name = {Species_name}\n')


AbEvol = AbEvol[:, :, :, :, iElm, :]
print('AbEvol.shape = ', AbEvol.shape)

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
nH_x = 10**3.0001
T_x = 10**5.60001
AbEvol_x = multidimensional_interpolate(rkpc, Lsh, nH, Temp, AbEvol, rkpc_x, Lsh_x, nH_x, T_x)

rkpc_x = 0.310001
AbEvol_x2= multidimensional_interpolate(rkpc, Lsh, nH, Temp, AbEvol, rkpc_x, Lsh_x, nH_x, T_x)

rkpc_x = 0.351
AbEvol_x3= multidimensional_interpolate(rkpc, Lsh, nH, Temp, AbEvol, rkpc_x, Lsh_x, nH_x, T_x)


TT = AbEvol[1, 3, 35, 13, :] # (N_rkpc, N_Lsh, N_nH, N_T, N_elm, N_time)
TT2 = AbEvol[2, 3, 35, 13, :] # (N_rkpc, N_Lsh, N_nH, N_T, N_elm, N_time)

plt.plot(t, np.log10(TT), color = 'k', linewidth = 3)
plt.plot(t, np.log10(TT2), color = 'grey', linewidth = 3)
plt.scatter(t, Abchim, s = 5, color = 'r')
plt.scatter(t, np.log10(AbEvol_x), s = 5, color = 'orange')
plt.scatter(t, np.log10(AbEvol_x2), s = 5, color = 'lime')
plt.scatter(t, np.log10(AbEvol_x3), s = 5, color = 'blue')

plt.savefig('fig2.png', bbox_inches = 'tight', dpi = 150)

plt.show()



