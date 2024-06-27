#import h5py
import numpy as np
import matplotlib.pyplot as plt
#from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp
import pickle
from data import *



#----- dT_dt
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
  
#----- grd_f_T
def grd_f_T(j, nH0, nHp, nHe0, nHep, nHepp, T):

  dx = np.zeros(6)
  delta = 1e-6
  dx[j] = delta
  grdx = (dT_dt(nH0+dx[0], nHp+dx[1], nHe0+dx[2], nHep+dx[3], nHepp+dx[4], T+dx[5]) - (dT_dt(nH0, nHp, nHe0, nHep, nHepp, T))) / delta

  return grdx




#----- dnH0_dt
def dnH0_dt(nH0, nHp, nHe0, nHep, nHepp, T):

  Tx = np.log10(T)
  ne = nHp + nHep + 2.0 * nHepp

  return 10**k2(Tx) * nHp * ne - 10**k1(Tx) * nH0 * ne

#----- grd_f_nH0
def grd_f_nH0(j, nH0, nHp, nHe0, nHep, nHepp, T):

  dx = np.zeros(6)
  delta = 1e-10
  dx[j] = delta
  grdx = (dnH0_dt(nH0+dx[0], nHp+dx[1], nHe0+dx[2], nHep+dx[3], nHepp+dx[4], T+dx[5]) - dnH0_dt(nH0, nHp, nHe0, nHep, nHepp, T)) / delta
  
  if j not in [0, 1]:
    return 0.0
  
  return grdx




#----- dnHp_dt
def dnHp_dt(nH0, nHp, nHe0, nHep, nHepp, T):
  
  Tx = np.log10(T)
  ne = nHp + nHep + 2.0 * nHepp
  
  return 10**k1(Tx) * nH0 * ne - 10**k2(Tx) * nHp * ne

#----- grd_f_nHp
def grd_f_nHp(j, nH0, nHp, nHe0, nHep, nHepp, T):

  dx = np.zeros(6)
  delta = 1e-6
  dx[j] = delta
  grdx = (dnHp_dt(nH0+dx[0], nHp+dx[1], nHe0+dx[2], nHep+dx[3], nHepp+dx[4], T+dx[5]) - dnHp_dt(nH0, nHp, nHe0, nHep, nHepp, T)) / delta
  
  if j not in [0, 1]:
    return 0.0
  
  return grdx




#----- dnHe0_dt
def dnHe0_dt(nH0, nHp, nHe0, nHep, nHepp, T):
  
  Tx = np.log10(T)
  ne = nHp + nHep + 2.0 * nHepp
  
  return 10**k5(Tx) * nHep * ne - 10**k3(Tx) * nHe0 * ne

#----- grd_f_nHe0
def grd_f_nHe0(j, nH0, nHp, nHe0, nHep, nHepp, T):

  dx = np.zeros(6)
  delta = 1e-6
  dx[j] = delta
  grdx = (dnHe0_dt(nH0+dx[0], nHp+dx[1], nHe0+dx[2], nHep+dx[3], nHepp+dx[4], T+dx[5]) - dnHe0_dt(nH0, nHp, nHe0, nHep, nHepp, T)) / delta
  
  if j not in [2, 3]:
    return 0.0
  
  return grdx



#----- dnHep_dt
def dnHep_dt(nH0, nHp, nHe0, nHep, nHepp, T):
  
  Tx = np.log10(T)
  ne = nHp + nHep + 2.0 * nHepp
  
  return 10**k6(Tx) * nHepp * ne + 10**k3(Tx) * nHe0 * ne - 10**k4(Tx) * nHep * ne - 10**k5(Tx) * nHep * ne

#----- grd_f_nHep
def grd_f_nHep(j, nH0, nHp, nHe0, nHep, nHepp, T):

  dx = np.zeros(6)
  delta = 1e-6
  dx[j] = delta
  grdx = (dnHep_dt(nH0+dx[0], nHp+dx[1], nHe0+dx[2], nHep+dx[3], nHepp+dx[4], T+dx[5]) - dnHep_dt(nH0, nHp, nHe0, nHep, nHepp, T)) / delta
  
  if j not in [2, 3, 4]:
    return 0.0
  
  return grdx




#----- dnHepp_dt
def dnHepp_dt(nH0, nHp, nHe0, nHep, nHepp, T):
  
  Tx = np.log10(T)
  ne = nHp + nHep + 2.0 * nHepp
  
  return 10**k4(Tx) * nHep * ne - 10**k6(Tx) * nHepp * ne
 
#----- grd_f_nHepp
def grd_f_nHepp(j, nH0, nHp, nHe0, nHep, nHepp, T):

  dx = np.zeros(6)
  delta = 1e-6
  dx[j] = delta
  grdx = (dnHepp_dt(nH0+dx[0], nHp+dx[1], nHe0+dx[2], nHep+dx[3], nHepp+dx[4], T+dx[5]) - dnHepp_dt(nH0, nHp, nHe0, nHep, nHepp, T)) / delta
  
  if j not in [3, 4]:
    return 0.0
  
  return grdx




nH0, nHp, nHe0, nHep, nHepp, T = 0.1, 0.9, 0.1, 0.3, 0.6, 20000



#def jacobian(nH0, nHp, nHe0, nHep, nHepp, T):

lst = []
lstX = []
for j in range(6):
  lst.append(grd_f_nH0(j, nH0, nHp, nHe0, nHep, nHepp, T))
lstX.append(lst)

lst = []
for j in range(6):
  lst.append(grd_f_nHp(j, nH0, nHp, nHe0, nHep, nHepp, T))
lstX.append(lst)

lst = []
for j in range(6):
  lst.append(grd_f_nHe0(j, nH0, nHp, nHe0, nHep, nHepp, T))
lstX.append(lst)

lst = []
for j in range(6):
  lst.append(grd_f_nHep(j, nH0, nHp, nHe0, nHep, nHepp, T))
lstX.append(lst)

lst = []
for j in range(6):
  lst.append(grd_f_nHepp(j, nH0, nHp, nHe0, nHep, nHepp, T))
lstX.append(lst)

lst = []
for j in range(6):
  lst.append(grd_f_T(j, nH0, nHp, nHe0, nHep, nHepp, T))
lstX.append(lst)


lstX = np.array(lstX)
print(lstX)
print('shape = ', lstX.shape)






#j = 5 # nH0 -------- 1-->nHp    2-->nHe0     3-->nHep    4-->nHepp     5-->T
#grdfT = grd_f_T(j, nH0, nHp, nHe0, nHep, nHepp, T)
#print(grdfT)

#j = 2 # nH0 -------- 1-->nHp    2-->nHe0     3-->nHep    4-->nHepp     5-->T
#grdfnHep = grd_f_nHep(j, nH0, nHp, nHe0, nHep, nHepp, T)
#print(grdfnH0)







