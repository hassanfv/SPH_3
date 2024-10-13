
import numpy as np
import pandas as pd


df = pd.read_csv('xOut.csv')

print(df.columns)
#['nH', 'rkpc', 'Lsh', 't10', 't9', 't8', 't7', 't6', 't5', 't4', 't3', 'min_T', 'TR1', 'TR2']

min_T = df['min_T']
TR1 = df['TR1']
TR2 = df['TR2']

nx = np.where(TR2 > 1e4)[0]


dfT = df.iloc[nx, :]
print(dfT.sort_values(by = 'nH'))




