
# Ref: Grassi et al - 2011

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.integrate import odeint


#===== LambColIonH0 ---> erg.s^-1.cm^-3
def LambColIonH0(T, ne, nH0):
  return 1.27e-21 * T**0.5 * np.exp(-157809.1/T) / (1.0 + (T/1e5)**0.5) * ne * nH0



#===== a_Hp
def a_Hp(T): 
  return 8.4e-11 / T**0.5 / (T/1e3)**0.2 / (1.0 + (T/1e6)**0.7)


#----------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
#------------------------------------- Reaction Rates (cm^3.s^-1) -----------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
# Reaction: (H0 + e ---> Hp + 2e) ::: H0 Collisional ionization
def k1(T):

  T = T * 8.61732814974493E-05 # convert Kelvin to eV.
  ln_T = np.log(T)
  
  k1_val = np.exp(-32.71396786 \
                  + 13.5365560 * ln_T \
                  - 5.73932875 * ln_T**2 \
                  + 1.56315498 * ln_T**3 \
                  - 0.287705600 * ln_T**4 \
                  + 3.48255977 * 10**-2 * ln_T**5 \
                  - 2.63197617 * 10**-3 * ln_T**6 \
                  + 1.11954395 * 10**-4 * ln_T**7 \
                  - 2.03914985 * 10**-6 * ln_T**8) 
  return k1_val
#--------------

# Reaction: (Hp + e ---> H0 + γ)  ::: photo-recombination of Hp and e.
def k2(T):

  k2_val = 2.753e-14 * (315614 / T)**1.5 * (1 + (115188 / T)**0.407)**-2.242 # Grassi made a mistake and wrote 1115188 instead of 115188!!

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

  g2_val = 13.5984 * 1.60218e-12 * k2(T)

  return g2_val


# Cooling due to H0 collisional excitation via collision with electrons ::: in erg.s^-1.cm^3 ::: will be multiplied by nH0 and ne later in the code.
def g3(T):

  g3_val = 7.50e-19 * np.exp(-118348/T) / (1.0 + (T/1e5)**0.5)
  
  return g3_val



T = 1e4
print(g3(T), g1(T), g2(T))
#s()


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
  
  print(np.log10(T), np.log10(Lam))
  #print(nH0, T)
  
  ntot = n_tot(nH, nH0)

  dT_dt = -1.0 * (gamma - 1.0) / kB / ntot * Lam
  
  return [dnH0_dt, dT_dt]



gamma = 5./3.
kB = 1.3807e-16

nH = 1.0

y0 = [1e-3, 1e6]

t = np.linspace(1, 120000, 10000) * 3.16e7
t_yrs = t/3.16e7

sol = odeint(func, y0 = y0, t = t, tfirst = True)

nH0 = sol[:, 0]
T = sol[:, 1]

print(T)



nHp = nH - nH0
ne = nHp

plt.subplot(1, 2, 1)
plt.scatter(t_yrs, (T), s = 5, color = 'k')
plt.ylim(3, 8)

plt.subplot(1, 2, 2)
plt.scatter(t_yrs, nHp, s = 5, color = 'k', label = 'nHp')
plt.scatter(t_yrs, nH0, s = 5, color = 'b', label = 'nH0')
plt.scatter(t_yrs, ne, s = 5, color = 'orange', label = 'ne')
plt.yscale('log')
plt.legend()

plt.savefig('myOnlyH.png')

plt.show()




