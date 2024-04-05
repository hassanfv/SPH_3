
import numpy as np
from photolibs3 import *
import matplotlib.pyplot as plt


#----- get_mu ----> ref: Townsend - 2009
def get_mu(amu, X, Z):
  return amu / (2.0 * X + 3.0 * (1.0 - X - Z)/4.0 + Z/2.0)


#----- get_tcool
def get_tcool(T_n, nH):
  Gam, Lamb_at_T_n = coolingHeatingRates(T_n, nH)
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

nHcgs = 1.0

#---- Test case ----
T = 50000
t_cool = get_tcool(T, nHcgs)
print(f't_cool = {t_cool:.3E} seconds')
#-------------------

Tgrid = np.logspace(1.0, 10.0, 1000)

print(Tgrid)

res = []

for T in Tgrid:

  Gam, Lam = coolingHeatingRates(T, nHcgs)

  res.append([Gam, Lam])

res = np.array(res)

Gam = res[:, 0]
Lam = res[:, 1]

plt.scatter(np.log10(Tgrid), np.log10(Lam), s = 1)
plt.show()





