
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


#--- INPUT ----
nHx = 1.0 # in LOG
rx = 0.1 # kpc
NHx = 19.0 # in LOG
Tx = 5.5

tol = 0.01
ndx_nH = (lognH >= nHx - tol) & (lognH <= nHx + tol)
ndx_r = (rkpc >= rx - tol) & (rkpc <= rx + tol)
ndx_NH = (logNH >= NHx - tol) & (logNH <= NHx + tol)

nx = np.where(ndx_nH & ndx_r & ndx_NH)[0]

H_tmp = Heat[nx]
C_tmp = Cool[nx]
T_tmp = logT[nx]

N_T = len(T_tmp) # as a unit test we can define N_T from the inputList code and then compare N_T and len(T_tmp), if they don't match raise error!!!!!

#---- Find the T smaller and bigger than Tx
for i in range(N_T):
  if Tx <= T_tmp[i]:
    nT = i
    break

if nT > 0:
  nT = nT - 1 # to make nT as the lower bound. We add 1 to nT later to get the upper bound.... we do this just for convenience!

print(T_tmp[nT], Tx, T_tmp[nT+1])

Tgrid = np.arange(2, 11.1, 0.1)

print(len(nx))
print(len(Tgrid))

print(f'Heat = {H_tmp[nT]},  Cool = {C_tmp[nT]}')
print(f'Heat - Cool = {10**H_tmp[nT] - 10**C_tmp[nT]}')


#-------------- From Gnat & Sternberg - 2007 ---------------------
df = pd.read_csv('table4a.dat', delim_whitespace=True, header=None)
T = df.iloc[:, 0].values
Lam = df.iloc[:, 3].values
#------------------------------------

coolx = np.log10(10**Cool[nx] / (10**nHx)**2)
heatx = np.log10(10**Heat[nx] / (10**nHx)**2)

plt.scatter(Tgrid, coolx, s = 1, color = 'blue')
plt.scatter(Tgrid, heatx, s = 1, color = 'red')

plt.scatter(T_tmp[nT], H_tmp[nT], s = 30, color = 'blue')
plt.plot(np.log10(T), np.log10(Lam), color = 'k')

plt.show()





