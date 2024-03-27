
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


df = pd.read_csv('HCmu.csv')
#['lognH', 'rkpc', 'logNHtot', 'logT', 'logHeating', 'logCooling', 'mu']

lognH = df['lognH'].values
rkpc = df['rkpc'].values
logNH = df['logNHtot'].values
logT = df['logT'].values

Heat = df['logHeating'].values
Cool = df['logCooling'].values


nHx = 4.0 # in LOG
rx = 0.3 # kpc
NHx = 20.0 # in LOG

tol = 0.01
ndx_nH = (lognH >= nHx - tol) & (lognH <= nHx + tol)
ndx_r = (rkpc >= rx - tol) & (rkpc <= rx + tol)
ndx_NH = (logNH >= NHx - tol) & (logNH <= NHx + tol)

nx = np.where(ndx_nH & ndx_r & ndx_NH)[0]

Tgrid = np.arange(2, 11.1, 0.1)

print(len(nx))
print(len(Tgrid))

plt.scatter(Tgrid, Cool[nx], s = 1, color = 'blue')
plt.scatter(Tgrid, Heat[nx], s = 1, color = 'red')

plt.show()


