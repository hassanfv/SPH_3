
import numpy as np
from photolibs3 import *
import matplotlib.pyplot as plt
import pandas as pd

# In this version we use cooling functions from Gnat & Sterberg - 2007 as my photolibs cooling did not consider metal cooling!!

#-------------- From Gnat & Sternberg - 2007 ---------------------
df = pd.read_csv('table4a.dat', delim_whitespace=True, header=None)
TGnat = df.iloc[:, 0].values
LGnat = df.iloc[:, 3].values

TGnat = np.log10(TGnat)
LGnat = np.log10(LGnat)
#------------------------------------


#----- get_tcool
def get_tcool2(Tx, nH):

  Tx = np.log10(Tx)

  nt = -1
  for i in range(len(TGnat)):
    if Tx < TGnat[i]:
      nt = i
      break

  if nt == -1:
    nt = len(TGnat) - 1

  if nt > 1:
    nt -= 1

  TL = TGnat[nt]
  TU = TGnat[nt+1]

  LL = LGnat[nt]
  LU = LGnat[nt+1]
    
  slope = (LU - LL) / (TU - TL)

  Lres = LL + slope * (Tx - TL)

  Lamb_at_Tx = 10**Lres

  rho = nH / X
  
  return kB * muelec * muH * 10**Tx / (gamma - 1.0) / rho / mu / Lamb_at_Tx


#----- get_mu ----> ref: Townsend - 2009
def get_mu(amu, X, Z):
  return amu / (2.0 * X + 3.0 * (1.0 - X - Z)/4.0 + Z/2.0)


#----- get_tcool
def get_tcool(T_n, nH):
  Gam, Lamb_at_T_n = coolingHeatingRates(T_n, nH) # This does not contain metal cooling which was considered in Townsend - 2009 paper!
  rho = nH / X
  return kB * muelec * muH * T_n / (gamma - 1.0) / rho / mu / Lamb_at_T_n


#----- update_T
def update_T(T_n, dt, tcool):
  return T_n * (1.0 - dt / tcool)



kB = 1.3807e-16  # cm2 g s-2 K-1
mH = 1.6726e-24 # gram
gamma = 5.0/3.0


amu = 1.0
X = 0.7
Z = 0.02

mu = get_mu(amu, X, Z)
print(mu)


nH = 100.0

muH = amu / X
muelec = 2.0 * amu / (1.0 + X)
print(mu, muH, muelec)

gJH0 = gJHe0 = gJHep = 0.0 # No radiation field!

nHcgs = 1000.0

#---- Test case ----
T = 1e6
t_cool = get_tcool2(T, nHcgs)
print(f't_cool = {t_cool:.3E} seconds')
print(f't_cool = {(t_cool/365.25/24./3600.):.3f} years')
#-------------------

dt_yrs = 2 # time-step
dt = dt_yrs * 365.25 * 24.0 * 3600.0

N = 2 * np.int32(t_cool / dt)

print('N = ', N)

res = []

t = 0.0

for i in range(N):
  
  tcool = get_tcool2(T, nHcgs)
  T = update_T(T, dt, tcool)
  
  if T < 1e4:
    T = 1e4
  
  t += dt
  
  res.append([t, T])

res = np.array(res)

t = res[:, 0]
T = res[:, 1]

plt.scatter(t/t_cool, T, s = 1)
plt.xscale('log')
plt.yscale('log')
plt.xlim(0.1, 2.0)

plt.show()





