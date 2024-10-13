
import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import time
import glob

TTA = time.time()

dirX = '/home/pc/Desktop/N_body_2024/SPH_3/cooling29June2024/Table_preparation_3/On_EC2_Oct12/pklFiles01/pklFiles/'

filez = glob.glob(dirX + '*.pkl')

Xres = []

for j in range(len(filez)):

  nam = filez[j] # Good for testing ---> nH_2.4_rkpc_0.81_Lsh_20.739.pkl

  print(j, nam)

  with open(nam, 'rb') as f:
    data = pickle.load(f)
    # ['TempEvol', 'AbundEvol', 'nH', 'rkpc', 'Lsh', 'Species_id', 'Species_name', 't_in_sec', 'nH_p', 'rkpc_p', 'Lsh_p']

  #print(data.keys())


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

  Xres.append(
              [nH_p, rkpc_p, Lsh_p,
              float(Res[0]), float(Res[1]), float(Res[2]), float(Res[3]),
              float(Res[4]), float(Res[5]), float(Res[6]), float(Res[7]),
              float(min(TEvol)), float(np.median(TEvol[-2000:-1000])),
              float(np.median(TEvol[-1000:]))]
            )

colz = ['nH', 'rkpc', 'Lsh', 't10', 't9', 't8', 't7', 't6', 't5', 't4', 't3', 'min_T', 'TR1', 'TR2']
df = pd.DataFrame(Xres)
df.columns = colz

print()
print(df)
print()

df = df.round(3)
df = df.map(lambda x: f"{x:.4E}" if x > 100 else f"{x:.3f}")
df.to_csv('xOut.csv', index = False)

print(f'Elapsed time = {time.time() - TTA:.2f} seconds.')




