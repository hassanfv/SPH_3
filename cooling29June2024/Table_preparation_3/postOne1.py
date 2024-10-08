
import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import time
import glob


filez = glob.glob('*.pkl')

j = 0

nam = filez[j] # Good for testing ---> nH_2.4_rkpc_0.81_Lsh_20.739.pkl

print(nam)

with open(nam, 'rb') as f:
  data = pickle.load(f)
  # ['TempEvol', 'AbundEvol', 'nH', 'rkpc', 'Lsh', 'Species_id', 'Species_name', 't_in_sec']

print(data.keys())

TEvol = data['TempEvol']
t_Arr_in_yrs = data['t_in_sec'] / 3600. / 24. / 365.25

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

#----- Plot Section --------
plt.scatter(t_Arr_in_yrs, np.log10(TEvol), s = 1)
for timex in Res:
  if timex > 0.0:
    plt.axvline(x = timex, linestyle = ':', color = 'b')
plt.show()
#--------------------------




