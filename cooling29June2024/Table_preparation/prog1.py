
import numpy as np
import h5py
import matplotlib.pyplot as plt


def find_nearest_indices(array, value):
    idx = np.searchsorted(array, value)
    if idx == 0:
        return 0, 1
    elif idx == len(array):
        return len(array) - 2, len(array) - 1
    else:
        return idx - 1, idx




rkpc = [0.6, 0.8, 1.0]
NH = [18, 19, 20]

nH = np.arange(-2., 2.01, 0.1) #!!!!!!!!!!!!!!!!!!!!!!! CHECK len to be consistent with CHIMES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



Res = np.zeros((len(rkpc), len(NH), 41, 5001))  # 41 is len(nH) and we can get it from TEvol.shape! 5001 is len(time)

print(Res.shape)

for i, rkpci in enumerate(rkpc):
  for j, NHj in enumerate(NH):
    nam = f'./hdf5_files/grid_noneq_evolution_NeuralNet_rkpc_{rkpci}_NH_{NHj}.hdf5'
    f = h5py.File(nam, 'r')
    
    TEvol = f['TemperatureEvolution'][:] # (1, 41, 1, 5001) ---> (T, nH, Z, t)
    
    t_Arr_in_sec = f['TimeArray_seconds'][:]
    t_Arr_in_yrs = t_Arr_in_sec / (3600. * 24. * 365.25)
    nH_arr = f['TableBins/Densities'][:]
    
    for k, nHk in enumerate(nH_arr):
      Res[i, j, k, :] = np.log10(TEvol[0, k, 0, :])

dictx = {'rkpc': rkpc, 'NH': NH, 'nH': nH_arr, 'T': Res}

#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------

# T ---> (3, 3, 41, 5001) ---> (rkpc, NH, nH, t)

rkpc_i = 0.75 # This is in kpc !
NH_i = 19.5 # This is in log !
nH_i = 1.85 # This is in log !

rkpc_idx_low, rkpc_idx_high = find_nearest_indices(rkpc, rkpc_i)
NH_idx_low, NH_idx_high = find_nearest_indices(NH, NH_i)
nH_idx_low, nH_idx_high = find_nearest_indices(nH, nH_i)

print(rkpc_idx_low, rkpc_idx_high)
print(NH_idx_low, NH_idx_high)
print(nH_idx_low, nH_idx_high)

Tarr = Res.copy()

T000 = Tarr[rkpc_idx_low, NH_idx_low, nH_idx_low]
T001 = Tarr[rkpc_idx_low, NH_idx_low, nH_idx_high]
T010 = Tarr[rkpc_idx_low, NH_idx_high, nH_idx_low]
T011 = Tarr[rkpc_idx_low, NH_idx_high, nH_idx_high]
T100 = Tarr[rkpc_idx_high, NH_idx_low, nH_idx_low]
T101 = Tarr[rkpc_idx_high, NH_idx_low, nH_idx_high]
T110 = Tarr[rkpc_idx_high, NH_idx_high, nH_idx_low]
T111 = Tarr[rkpc_idx_high, NH_idx_high, nH_idx_high]


xd = (rkpc_i - rkpc[rkpc_idx_low]) / (rkpc[rkpc_idx_high] - rkpc[rkpc_idx_low])
yd = (NH_i - NH[NH_idx_low]) / (NH[NH_idx_high] - NH[NH_idx_low])
zd = (nH_i - nH[nH_idx_low]) / (nH[nH_idx_high] - nH[nH_idx_low])

T_interp = (
            T000 * (1 - xd) * (1 - yd) * (1 - zd) +
            T001 * (1 - xd) * (1 - yd) * zd +
            T010 * (1 - xd) * yd * (1 - zd) +
            T011 * (1 - xd) * yd * zd +
            T100 * xd * (1 - yd) * (1 - zd) +
            T101 * xd * (1 - yd) * zd +
            T110 * xd * yd * (1 - zd) +
            T111 * xd * yd * zd
           )

print()
print(T_interp)

print()
print(t_Arr_in_yrs)

plt.scatter(t_Arr_in_yrs, Tarr[rkpc_idx_low, NH_idx_low, nH_idx_low, :], color = 'b', s = 5)
plt.scatter(t_Arr_in_yrs, T_interp, color = 'r', s = 5)
plt.scatter(t_Arr_in_yrs, Tarr[rkpc_idx_high, NH_idx_high, nH_idx_high, :], color = 'green', s = 5)

plt.title(f'{rkpc_i}')

plt.ylim(3.9, 5.1)

#plt.yscale('log')

plt.show()


