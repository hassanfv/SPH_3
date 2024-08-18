
import numpy as np
import h5py
import matplotlib.pyplot as plt
import pickle


def find_nearest_indices(values, target):
    # Find indices of the nearest lower and higher values
    idx_below = np.max(np.where(values <= target)[0])
    idx_above = np.min(np.where(values >= target)[0])
    return idx_below, idx_above

def linear_interpolate(x, y, x0):
    # Linear interpolation between two points (x0 lies between x1 and x2)
    x1, x2 = x
    y1, y2 = y
    return y1 + (y2 - y1) * (x0 - x1) / (x2 - x1)

def bilinear_interpolate(rkpc, NH, Tarr, rkpc_i, NH_i, time_index):
    # Find nearest indices in both dimensions
    rkpc_idx1, rkpc_idx2 = find_nearest_indices(rkpc, rkpc_i)
    NH_idx1, NH_idx2 = find_nearest_indices(NH, NH_i)
    
    # Interpolate along rkpc for both NH values at specific time index
    temp_at_NH1 = linear_interpolate(
        [rkpc[rkpc_idx1], rkpc[rkpc_idx2]],
        [Tarr[rkpc_idx1, NH_idx1, time_index], Tarr[rkpc_idx2, NH_idx1, time_index]],
        rkpc_i
    )
    temp_at_NH2 = linear_interpolate(
        [rkpc[rkpc_idx1], rkpc[rkpc_idx2]],
        [Tarr[rkpc_idx1, NH_idx2, time_index], Tarr[rkpc_idx2, NH_idx2, time_index]],
        rkpc_i
    )
    
    # Interpolate between the results along NH
    interpolated_temp = linear_interpolate(
        [NH[NH_idx1], NH[NH_idx2]],
        [temp_at_NH1, temp_at_NH2],
        NH_i
    )
    return interpolated_temp


def generate_temperature_evolution(rkpc, NH, Tarr, rkpc_i, NH_i, total_time):
    return np.array([bilinear_interpolate(rkpc, NH, Tarr, rkpc_i, NH_i, t) for t in range(total_time)])






rkpc = np.array([0.50, 0.55, 0.60])
NH = np.array([19.0, 19.2, 19.4])

nH = np.arange(-2., 2.01, 0.1) #!!!!!!!!!!!!!!!!!!!!!!! CHECK len to be consistent with CHIMES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Res = np.zeros((len(rkpc), len(NH), 5001))  # 41 is len(nH) and we can get it from TEvol.shape! 5001 is len(time)

print(Res.shape)

for i, rkpci in enumerate(rkpc):
  for j, NHj in enumerate(NH):
    nam = f'./hdf5_files/grid_noneq_evolution_NeuralNet_rkpc_{rkpci:.2f}_NH_{NHj:.1f}.hdf5'
    f = h5py.File(nam, 'r')
    
    TEvol = f['TemperatureEvolution'][:] # (1, 41, 1, 5001) ---> (T, nH, Z, t)
    
    t_Arr_in_sec = f['TimeArray_seconds'][:]
    t_Arr_in_yrs = t_Arr_in_sec / (3600. * 24. * 365.25)
    nH_arr = f['TableBins/Densities'][:]

    Res[i, j, :] = np.log10(TEvol[0, -1, 0, :]) # -1 means we take the last nH


#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------

# T ---> (3, 3, 41, 5001) ---> (rkpc, NH, nH, t)

rkpc_i = 0.58 # This is in kpc !
NH_i = 19.35 # This is in log !

Tarr = Res.copy()

total_time = 5001
T_interp = generate_temperature_evolution(rkpc, NH, Tarr, rkpc_i, NH_i, total_time)



print()
print('T_interp = ', T_interp)

plt.scatter(t_Arr_in_yrs, T_interp, color = 'r', s = 5)

plt.title(f'{rkpc_i}')

plt.ylim(3.9, 5.1)
#plt.yscale('log')

plt.show()




