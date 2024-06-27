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



#----- jacobian
def jacobian(nH0, nHp, nHe0, nHep, nHepp, T):

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

  return lstX


# Total cooling rate in erg.s^-1.cm^-3 ===>NOTE that Hep and Hepp is excluded in free-free here as we are only considering H in this code!
def Lambda(T, nH0, nHp, nHe0, nHep, nHepp):

  Tx = np.log10(T)
  
  ne = nHp + nHep + 2.0 * nHepp
  
  Lamb = (
           10**g1(Tx) * ne * nH0  # H0
         + 10**g2(Tx)  * ne * nHp # Hp
         + 10**g3(Tx) * nHe0 * ne # He0 
         + 10**g4(Tx) * nHep * ne # Hep 
         + 10**g5(Tx) * nHepp * ne# Hepp
        )
  
  return Lamb


#-----------------------------------------
#------ Solution of the ODE Section ------
#-----------------------------------------
def func(t, y):

  nH0, nHp, nHe0, nHep, nHepp, T = y
  
  Tx= np.log10(T)
  
  ne = nHp + nHep + 2.0 * nHepp
  ntot = nH0 + nHp + nHe0 + nHep + nHepp + ne
  
  dnH0_dt = 10**k2(Tx) * nHp * ne - 10**k1(Tx) * nH0 * ne
  
  dnH0_dt = 10**k1(Tx) * nH0 * ne - 10**k2(Tx) * nHp * ne
  
  dnHe0_dt = 10**k5(Tx) * nHep * ne - 10**k3(Tx) * nHe0 * ne
  
  dnHep_dt = 10**k6(Tx) * nHepp * ne + 10**k3(Tx) * nHe0 * ne - 10**k4(Tx) * nHep * ne - 10**k5(Tx) * nHep * ne
  
  dnHep_dt = 10**k4(Tx) * nHep * ne - 10**k6(Tx) * nHepp * ne
  
  Lamb = Lambda(T, nH0, nHp, nHe0, nHep, nHepp)

  dT_dt = -1.0 * (gamma - 1.0) / kB / ntot * Lamb
  
  return [dnH0_dt, dnHp_dt, dnHe0_dt, dnHep_dt, dnHepp_dt, dT_dt]





nH = 1000.0

He_solar = 10**(-1.07)
nHe = He_solar * nH
print('nHe (cm^-3) = ', nHe)

nH0_i = 1.0
nHp_i = 999.0

nHe0_i = 0.01 * nHe
nHep_i = 0.1 * nHe
nHepp_i = nHe - nHe0_i - nHep_i

y0 = [nH0_i, nHp_i, nHe0_i, nHep_i, nHepp_i, 1e6]

t_span = (1*3.16e7, 20000*3.16e7)

solution = solve_ivp(func, t_span, y0, method='LSODA', jac=jacobian, dense_output=True)

t = np.linspace(t_span[0], t_span[1], 10000) # This 10000 is not years, it is the number of points in linspace !!!!
y = solution.sol(t)

print(y.shape)










