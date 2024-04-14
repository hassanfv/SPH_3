
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


#--------------------------------
file_path = 'CarbonCoolingRates.txt'
df = pd.read_csv(file_path, delim_whitespace=True, skiprows=19) # Adjust skiprows as necessary
df.columns = ['Temp', 'C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'CCIE-only-C']
Temp = df['Temp'].values
C0_LAMBDA = df['C0'].values
C1_LAMBDA = df['C1'].values
C2_LAMBDA = df['C2'].values
C3_LAMBDA = df['C3'].values
C4_LAMBDA = df['C4'].values
C5_LAMBDA = df['C5'].values
C6_LAMBDA = df['C6'].values
#---------------------------------


T = 14569.0

#----- closestx
def closestx(arr, x0):
  
  nt = np.argmin(np.abs(arr - x0))
  
  return arr[nt]


print(closestx(Temp, T))



