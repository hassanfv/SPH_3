
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






rkpc = np.array([0.50, 0.60, 0.70])
NH = np.array([20.5, 21.0, 21.5])

#nH = np.arange(-2., 2.01, 0.1) #!!!!!!!!!!!!!!!!!!!!!!! CHECK len to be consistent with CHIMES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Res = np.zeros((len(rkpc), len(NH), 5001))  # 41 is len(nH) and we can get it from TEvol.shape! 5001 is len(time)

print(Res.shape)

for i, rkpci in enumerate(rkpc):
  for j, NHj in enumerate(NH):
    nam = f'./hdf5_files/grid_noneq_evolution_NeuralNetX_r_kpc_{rkpci:.2f}_NH_{NHj:.1f}.hdf5'
    f = h5py.File(nam, 'r')
    
    TEvol = f['TemperatureEvolution'][:] # (1, 41, 1, 5001) ---> (T, nH, Z, t)
    
    t_Arr_in_sec = f['TimeArray_seconds'][:]
    t_Arr_in_yrs = t_Arr_in_sec / (3600. * 24. * 365.25)
    #nH_arr = f['TableBins/Densities'][:]

    Res[i, j, :] = np.log10(TEvol[0, -1, 0, :]) # -1 means we take the last nH


#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------

rkpc_i = 0.55 # This is in kpc !
NH_i = 21.25 # This is in log !

Tarr = Res.copy()

total_time = 5001
T_interp = generate_temperature_evolution(rkpc, NH, Tarr, rkpc_i, NH_i, total_time)


#----- The file is created in test1.py --------- Check if tmp.pkl is updated for this current data in prog1.py !!!!
with open('tmp.pkl', 'rb') as f:
  tmp = pickle.load(f)
tt = tmp['t']
TT = tmp['T']
#-----------------------------------------------

plt.scatter(t_Arr_in_yrs, T_interp, color = 'r', s = 5)
plt.scatter(tt, TT, color = 'gold', s = 5) #--------- Check if tmp.pkl is updated for this currect data in prog1.py !!!!

plt.title(f'{rkpc_i}')

plt.ylim(3.75, 5.1)

plt.show()




