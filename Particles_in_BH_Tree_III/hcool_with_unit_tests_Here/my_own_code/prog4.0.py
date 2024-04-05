
# This version uses odeint but prog4.1.py version uses solve_ivp.
# This works well for only Hydrogen atom as the coolant. Free-free is included in prog5.0.py

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.integrate import odeint

#----------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
#------------------------------------- Reaction Rates (cm^3.s^-1) -----------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
# Reaction: (H0 + e ---> Hp + 2e) ::: H0 Collisional ionization
def k1(T):
  
  k1_val = 5.85e-11 * T**0.5 * np.exp(-157809.1/T) / (1.0 + (T/1e5)**0.5)

  return k1_val
#--------------

# Reaction: (Hp + e ---> H0 + γ)  ::: photo-recombination of Hp and e.
def k2(T):

  k2_val = 8.4e-11 / T**0.5 / (T/1e3)**0.2 / (1.0 + (T/1e6)**0.7)

  return k2_val
#--------------


#----------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
#------------------------------------- Cooling Rates (erg.s^-1.cm^3) --------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
# Cooling rate due to H collisional ionization: (H0 + e ---> Hp + 2e)  ::: in erg.s^-1.cm^3  ::: will be multiplied by nH0 and ne later in the code.
def g1(T):
  
  g1_val = 13.5984 * 1.60218e-12 * k1(T)
  
  return g1_val


# Cooling rate due to Hp photo-recombination: (Hp + e ---> H0 + γ)  ::: in erg.s^-1.cm^3  ::: will be multiplied by nHp and ne later in the code.
def g2(T):

  g2_val = 8.70e-27 * T**0.5 / (T/1e3)**0.2 / (1.0 + (T/1e6)**0.7)

  return g2_val


# Cooling due to H0 collisional excitation via collision with electrons ::: in erg.s^-1.cm^3 ::: will be multiplied by nH0 and ne later in the code.
def g3(T):

  g3_val = 7.50e-19 * np.exp(-118348/T) / (1.0 + (T/1e5)**0.5)
  
  return g3_val


# Total cooling rate in erg.s^-1.cm^-3
def Lambda(T, nH, nH0):

  nHp = nH - nH0
  ne = nHp
  
  Lamb = (g1(T) * ne * nH0 # collisional ionization of hydrogen ::::::: (H0 + e ---> Hp + 2e).
       + g2(T) * ne * nHp # photo-recombination of Hp with electron :: (Hp + e ---> H0 + γ).
       + g3(T) * ne * nH0) # collisional excitaion of H0.
  
  return Lamb




def n_tot(nH, nH0):
  
  nHp = nH - nH0
  ne = nHp
  
  ntot = nH0 + nHp + ne
  
  return ntot



#-----------------------------------------
#------ Solution of the ODE Section ------
#-----------------------------------------
def func(t, y):

  nH0, T = y
  
  nHp = nH - nH0
  ne = nHp
  
  dnH0_dt = k2(T) * nHp * ne - k1(T) * nH0 * ne

  Lam = Lambda(T, nH, nH0)
  
  ntot = n_tot(nH, nH0)

  dT_dt = -1.0 * (gamma - 1.0) / kB / ntot * Lam
  
  return [dnH0_dt, dT_dt]



gamma = 5./3.
kB = 1.3807e-16

nH = 100.0

y0 = [2e-4, 1e6]

t = np.linspace(1, 120000, 10000) * 3.16e7
t_yrs = t/3.16e7

sol = odeint(func, y0 = y0, t = t, tfirst = True)

nH0 = sol[:, 0]
T = sol[:, 1]

print('T = ', T)

nHp = nH - nH0
ne = nHp

plt.figure(figsize = (12, 6))

plt.subplot(1, 2, 1)
plt.scatter(t_yrs, np.log10(T), s = 5, color = 'k')
plt.ylim(3, 8)

plt.subplot(1, 2, 2)
plt.scatter(t_yrs, nHp, s = 5, color = 'k', label = 'nHp')
plt.scatter(t_yrs, nH0, s = 5, color = 'b', label = 'nH0')
plt.scatter(t_yrs, ne, s = 5, color = 'orange', label = 'ne')
plt.yscale('log')
plt.title('odeint')
plt.legend()

plt.savefig('myOnlyH.png')

plt.show()




