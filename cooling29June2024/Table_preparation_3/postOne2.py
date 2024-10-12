
import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import time
import glob


filez = glob.glob('./pklFiles/*.pkl')

j = 1

nam = filez[j] # Good for testing ---> nH_2.4_rkpc_0.81_Lsh_20.739.pkl

print(nam)

with open(nam, 'rb') as f:
  data = pickle.load(f)
  # ['TempEvol', 'AbundEvol', 'nH', 'rkpc', 'Lsh', 'Species_id', 'Species_name', 't_in_sec', 'nH_p', 'rkpc_p', 'Lsh_p']

print(data.keys())


TEvol = data['TempEvol']
t_Arr_in_yrs = data['t_in_sec'] / 3600. / 24. / 365.25

nH_p = data['nH_p']
rkpc_p = data['rkpc_p']
Lsh_p = data['Lsh_p']

Res = np.zeros(9)

checker1 = checker2 = checker3 = checker4 = checker5 = checker6 = checker7 = checker8 = 0

#=========================
for i, T in enumerate(TEvol):
  if (np.log10(T) < 10.0) and (checker1 == 0):
    Res[0] = t_Arr_in_yrs[i]
    checker1 = 1
  
  if (np.log10(T) < 9.0) and (checker2 == 0):
    Res[1] = t_Arr_in_yrs[i]
    checker2 = 1
    
  if (np.log10(T) < 8.0) and (checker3 == 0):
    Res[2] = t_Arr_in_yrs[i]
    checker3 = 1
  
  if (np.log10(T) < 7.0) and (checker4 == 0):
    Res[3] = t_Arr_in_yrs[i]
    checker4 = 1
    
  if (np.log10(T) < 6.0) and (checker5 == 0):
    Res[4] = t_Arr_in_yrs[i]
    checker5 = 1
  
  if (np.log10(T) < 5.0) and (checker6 == 0):
    Res[5] = t_Arr_in_yrs[i]
    checker6 = 1
    
  if (np.log10(T) < 4.0) and (checker7 == 0):
    Res[6] = t_Arr_in_yrs[i]
    checker7 = 1
  
  if (np.log10(T) < 3.0) and (checker8 == 0):
    Res[7] = t_Arr_in_yrs[i]
    checker8 = 1
#=========================

print(f'Years in yrs for different T = {Res} ')

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# Include (Res, min(T), T_R1, T_R2)

aa = [0.0, 0.0, 0.0]
#aa += [Res[0], Res[1], Res[2], Res[3], Res[4], Res[5], Res[6], Res[7], min(TEvol), np.median(TEvol[-2000:-1000]), np.median(TEvol[-1000:])]
aa += [
    float(Res[0]), float(Res[1]), float(Res[2]), float(Res[3]),
    float(Res[4]), float(Res[5]), float(Res[6]), float(Res[7]),
    float(min(TEvol)), float(np.median(TEvol[-2000:-1000])),
    float(np.median(TEvol[-1000:]))]

pd.options.display.float_format = lambda x: f'{x:.2e}' if abs(x) >= 1e6 else f'{x:.3f}'

colz = ['nH', 'rkpc', 'Lsh', 'T10', 'T9', 'T8', 'T7', 'T6', 'T5', 'T4', 'T3', 'min_T', 'TR1', 'TR2']

df = pd.DataFrame(aa)
df = df.T
df.columns = colz

print()
print(df)
print()


nHG = np.arange(-4.0, 4.01, 0.1)
rkpcG = np.arange(0.01, 1.02, 0.1)
LshG = np.arange(0.0, 2.51, 0.25)

nHstep = 0.1
rkpcstep = 0.1
Lshstep = 0.25

#nH_p = 0.4
#rkpc_p = 0.7
#Lsh_p = 0.75

ndx_nH = round((nH_p - nHG[0]) / nHstep)
ndx_rkpc = round((rkpc_p - rkpcG[0]) / rkpcstep)
ndx_Lsh = round((Lsh_p - LshG[0]) / Lshstep)

print(f'nH_p = {nH_p},  nHG[ndx_nH] = {nHG[ndx_nH]}')
print(f'rkpc_p = {rkpc_p},  rkpcG[ndx_rkpc] = {rkpcG[ndx_rkpc]}')
print(f'Lsh_p = {Lsh_p},  LshG[ndx_Lsh] = {LshG[ndx_Lsh]}')


#----- Plot Section --------
plt.scatter(t_Arr_in_yrs, np.log10(TEvol), s = 1)

for timex in Res:
  if timex > 0.0:
    plt.axvline(x = timex, linestyle = ':', color = 'b')

plt.axvline(x = t_Arr_in_yrs[-1000], linestyle = ':', color = 'red')
plt.axvline(x = t_Arr_in_yrs[-2000], linestyle = ':', color = 'red')

plt.show()
#--------------------------




