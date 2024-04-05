
import numpy as np



#----- get_mu ----> ref: Townsend - 2009
def get_mu(amu, X, Z):
  return amu / (2.0 * X + 3.0 * (1.0 - X - Z)/4.0 + Z/2.0)


def get_tcool(T_n):
  
  Gam_T_n = getGam(T_n, nH)
  return kB * muelec * muH * T_n / (gamma - 1.0) / rho / mu / Gam_T_n
  




amu = 1.0
X = 0.7
Z = 0.02

mu = get_mu(amu, X, Z)
print(mu)


nH = 100.0

muH = amu / X
muelec = 2.0 * amu / (1.0 + X)
print(mu, muH, muelec)






