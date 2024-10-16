
import numpy as np
import pandas as pd


df = pd.read_csv('xOut.csv')

print(df.columns)
#['nH', 'rkpc', 'Lsh', 't10', 't9', 't8', 't7', 't6', 't5', 't4', 't3', 'min_T', 'TR1', 'TR2']

#df['Lsh'] = 10**df['Lsh'] * 3.086e18

min_T = df['min_T']
TR1 = df['TR1']
TR2 = df['TR2']

#nx = np.where(TR2 <= 1e4)[0]
#nx = np.where((min_T > 9700.) & (min_T < 10000.))[0]
#nx = np.where((min_T > 10000.) & (TR2 < 11000.))[0]

nx = np.where(min_T <= 1e4)[0]

print((nx))


dfT = df.iloc[nx, :]

dfX = dfT.sort_values(by = ['nH', 'rkpc', 'Lsh'])

dfX.to_csv('useForMuRun.csv', index = False)

print(dfX)




