#import h5py
import numpy as np
import matplotlib.pyplot as plt
#from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp
import pickle
from data import *


def dT_dt(nH0, nHp, nHe0, nHep, nHepp, T):

  Tx = np.log10(T)
  
  ne = nHp + nHep + 2.0 * nHepp
  ntot = nH0 + nHp + nHe0 + nHep + nHepp + ne
  
  Lamb = (
           10**g1(Tx) * ne * nH0  # H0
         + 10**g2(Tx)  * ne * nHp # Hp
         + 10**g3(Tx) * nHe0 * ne # He0 
         + 10**g4(Tx) * nHep * ne # Hep 
         + 10**g5(Tx) * nHepp * ne# Hepp
        )
  
  dT_dt = -1.0 * (gamma - 1.0) / kB / ntot * Lamb
  
  return dT_dt


nH0, nHp, nHe0, nHep, nHepp, T = 0.1, 0.9, 0.1, 0.3, 0.6, 20000

ne = nHp + nHep + 2.0 * nHepp
ntot = nH0 + nHp + nHe0 + nHep + nHepp + ne
Tx = np.log10(T)
Lamb = (
           10**g1(Tx) * ne * nH0  # H0
         + 10**g2(Tx)  * ne * nHp # Hp
         + 10**g3(Tx) * nHe0 * ne # He0 
         + 10**g4(Tx) * nHep * ne # Hep 
         + 10**g5(Tx) * nHepp * ne# Hepp
        )


grd = -(gamma - 1.) / kB * ((-Lamb/ntot/ntot) + (10**g1(Tx) * ne / ntot)) # FOR nH0: direct derivation - i.e. using pen and paper!
print(grd)


def grd_f_T(j, nH0, nHp, nHe0, nHep, nHepp, T):

  dx = np.zeros(6)
  delta = 1e-6
  dx[j] = delta
  grdx = (dT_dt(nH0+dx[0], nHp+dx[1], nHe0+dx[2], nHep+dx[3], nHepp+dx[4], T+dx[5]) - (dT_dt(nH0, nHp, nHe0, nHep, nHepp, T))) / delta

  return grdx


j = 0 # nH0 -------- 1-->nHp    2-->nHe0     3-->nHep    4-->nHepp     5-->T
grdx = grd_f_T(j, nH0, nHp, nHe0, nHep, nHepp, T)

print(grdx)





