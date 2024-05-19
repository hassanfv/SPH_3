#import h5py
import numpy as np
import matplotlib.pyplot as plt
#from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp
import pickle
from data import *


def f_T(nH0, nHp, nHe0, nHep, nHepp, T):

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


grd = -(gamma - 1.) / kB * ((-Lamb/ntot/ntot) + (10**g1(Tx) * ne / ntot)) # direct derivation - i.e. using pen and paper!
print(grd)



delta = 1e-6
grdx = (f_T(nH0+delta, nHp, nHe0, nHep, nHepp, T) - (f_T(nH0, nHp, nHe0, nHep, nHepp, T))) / delta

print(grdx)


#fT = f_T(nH0, nHp, nHe0, nHep, nHepp, T)
#print(fT)

