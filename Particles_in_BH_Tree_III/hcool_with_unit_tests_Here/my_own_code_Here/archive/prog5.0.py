
# This version also includes Free-Free cooling (So we have cooling due to Hydrogen atoms and free-free emission).
# This version uses solve_ivp.

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


#!!! ::: will be multiplied by (nHp+nHep+nHepp) and ne later in the code.
def g4H(T):

  gff = 1.1 + 0.34 * np.exp(-(5.5 - np.log10(T))**2/3.0)
  g4H_val = 1.42e-27 * gff * T**0.5
  
  return g4H_val



# Total cooling rate in erg.s^-1.cm^-3 ===>NOTE that Hep and Hepp is excluded in free-free here as we are only considering H in this code!
def Lambda(T, nH, nH0):

  nHp = nH - nH0
  ne = nHp
  
  nHep = 0.0
  nHepp = 0.0
  
  Lamb = (g1(T) * ne * nH0 # collisional ionization of hydrogen ::::::: (H0 + e ---> Hp + 2e).
       + g2(T)  * ne * nHp # photo-recombination of Hp with electron :: (Hp + e ---> H0 + γ).
       + g3(T)  * ne * nH0 # collisional excitaion of H0.
       + g4H(T) * ne * (nHp + nHep + 4.0 * nHepp))
  
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

#t = np.linspace(1, 120000, 10000) * 3.16e7
#t_yrs = t/3.16e7

t_span = (1*3.16e7, 120000*3.16e7)

#sol = odeint(func, y0 = y0, t = t, tfirst = True)

solution = solve_ivp(func, t_span, y0, method='LSODA', dense_output=True)

t = np.linspace(t_span[0], t_span[1], 10000)
y = solution.sol(t)

print(y.shape)

t_yrs = t / 3.16e7

nH0 = y[0, :]
T = y[1, :]

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
plt.title('solve_ivp')
plt.legend()

plt.savefig('myOnlyH.png')

plt.show()



