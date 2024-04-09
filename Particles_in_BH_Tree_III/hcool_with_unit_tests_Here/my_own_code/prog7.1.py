
# In this version (i.e. prog7.1.py) we also include cooling rate plot!
# This version (i.e. prog7.0.py) uses kB * T * k2(T) as the cooling due to recombination of Hp. We can convert recomb. rate to cooling rate by multiplying kB*T!
# This version (i.e. prog6.0.py) also includs cooling due to He atom and ions.
# This version also includes Free-Free cooling (So we have cooling due to Hydrogen atoms and free-free emission).
# This version uses solve_ivp.

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.integrate import odeint
import pickle

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

# Reaction: (He0 + e ---> Hep + 2e) ::: He0 Collisional ionization 
def k3(T):

  k3_val = 2.38e-11 * T**0.5 * np.exp(-285335.4/T) / (1.0 + (T/1e5)**0.5)
  
  return k3_val
#--------------

# Reaction: (Hep + e ---> Hepp + 2e) ::: Hep Collisional ionization 
def k4(T):

  k4_val = 5.68e-12 * T**0.5 * np.exp(-631515./T) / (1.0 + (T/1e5)**0.5)
  
  return k4_val
#--------------

# Reaction: (Hep + e ---> He0 + γ)  ::: photo-recombination of Hep and e.
def k5(T):

  k5_val = 1.5e-10 / T**0.6353
  
  return k5_val
#--------------

# Reaction: (Hepp + e ---> Hep + γ)  ::: photo-recombination of Hepp and e.
def k6(T):

  k6_val = 3.36e-10 / T**0.5 / (T/1e3)**0.2 / (1. + (T/1e6)**0.7)
  
  return k6_val
#--------------

#Reaction: Hep di-electric recombination (Hep + e ---> He0 + γ)
def k7(T):

  k7_val = 1.9e-3 / T**1.5 * np.exp(-470000./T) * (1.0 + 0.3 * np.exp(-94000./T))
  
  return k7_val
#--------------





#----------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
#------------------------------------- Cooling Rates (erg.s^-1.cm^3) --------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
# Cooling rate due to H collisional ionization: (H0 + e ---> Hp + 2e)  ::: in erg.s^-1.cm^3  ::: will be multiplied by nH0 and ne later in the code.
def g1(T):

  g1_val = 1.27e-21 * T**0.5 * np.exp(-157809.1/T) / (1.0 + (T/1e5)**0.5)
  
  return g1_val


# Cooling rate due to Hp photo-recombination: (Hp + e ---> H0 + γ)  ::: in erg.s^-1.cm^3  ::: will be multiplied by nHp and ne later in the code.
def g2(T):

  #g2_val = 8.70e-27 * T**0.5 / (T/1e3)**0.2 / (1.0 + (T/1e6)**0.7)
  g2_val = kB * T * k2(T) # See LaMothe and Ferland - 2001, page 1, the lines below equation 1!

  return g2_val


# Cooling due to H0 collisional excitation via collision with electrons ::: in erg.s^-1.cm^3 ::: will be multiplied by nH0 and ne later in the code.
def g3(T):

  g3_val = 7.50e-19 * np.exp(-118348/T) / (1.0 + (T/1e5)**0.5)
  
  return g3_val


# Cooling due to free-free emission ::: will be multiplied by (nHp+nHep+nHepp)*ne later in the code.
def g4(T):

  gff = 1.1 + 0.34 * np.exp(-(5.5 - np.log10(T))**2/3.0)
  g4_val = 1.42e-27 * gff * T**0.5
  
  return g4_val


# Cooling due to Hep collisional excitation ::: will be multiplied by nHep and ne later in the code
def g5(T):

  g5_val = 5.54e-17 / T**0.397 * np.exp(-473638./T) / (1. + (T/1e5)**0.5)
  
  return g5_val


# Cooling due to He0 collisional ioniozation ::: will be multiplied by nHe0 and ne later in the code
def g6(T):
  
  g6_val = 9.38e-22 * T**0.5 * np.exp(-285335.4/T) / (1. + (T/1e5)**0.5)
  
  return g6_val


# Cooling via Hep collisional ionization ::: will be multiplied by nHep and ne later in the code
def g7(T):

  g7_val = 4.95e-22 * T**0.5 * np.exp(-631515./T) / (1. + (T/1e5)**0.5)
  
  return g7_val


# Cooling via Hep recombination with electron ::: will be multiplied by nHep and ne later in the code
def g8(T):

  g8_val = 1.55e-26 * T**0.3647
  
  return g8_val


# Cooling via Hepp recombination with electron ::: will be multiplied by nHepp and ne later in the code
def g9(T):
  
  g9_val = 3.48e-26 * T**0.5 / (T/1e3)**0.2 / (1. + (T/1e6)**0.7)
  
  return g9_val


# Cooling via di-electronic recombination of Hep with electrons ::: will be multiplied by nHep and ne later in the code
def g10(T):
  
  g10_val = 1.24e-13 / T**1.5 * np.exp(-470000./T) * (1. + 0.3 * np.exp(-94000./T))

  return g10_val





# Total cooling rate in erg.s^-1.cm^-3 ===>NOTE that Hep and Hepp is excluded in free-free here as we are only considering H in this code!
def Lambda(T, nH, nH0, nHe0, nHep):

  nHp = nH - nH0
  y = Y / 4.0 / (1.0 - Y) # Eq. 32 in Katz & Weinberg - 1996 (Y = 0.24, i.e. He mass fraction. X = 1.0 - Y)
  nHepp = y * nH - nHe0 - nHep
  
  ne = nHp + nHep + 2.0 * nHepp
  
  Lamb = (
           g1(T) * ne * nH0 # collisional ionization of hydrogen ::::::: (H0 + e ---> Hp + 2e).
         + g2(T)  * ne * nHp # photo-recombination of Hp with electron :: (Hp + e ---> H0 + γ).
         + g3(T)  * ne * nH0 # collisional excitaion of H0.
         + g4(T) * ne * (nHp + nHep + 4.0 * nHepp) # free-free emission
         + g5(T) * nHep * ne # collisional excitation of Hep.
         + g6(T) * nHe0 * ne # He0 collisional ionization
         + g7(T) * nHep * ne # Hep collisional ionization
         + g8(T) * nHep * ne # Hep recombination to He0
         + g9(T) * nHepp * ne# Hepp recombination to Hep
         + g10(T) * nHep * ne# Hep di-electric recombination to He0
        )
  
  return Lamb



#----- n_tot
def n_tot(nH, nH0, nHe0, nHep):
  
  nHp = nH - nH0

  y = Y / 4.0 / (1.0 - Y) # Eq. 32 in Katz & Weinberg - 1996 (Y = 0.24, i.e. He mass fraction. X = 1.0 - Y)
  nHepp = y * nH - nHe0 - nHep
  
  ne = nHp + nHep + 2.0 * nHepp
  
  ntot = nH0 + nHp + nHe0 + nHep + nHepp + ne
  
  return ntot



#-----------------------------------------
#------ Solution of the ODE Section ------
#-----------------------------------------
def func(t, y):

  nH0, nHe0, nHep, T = y
  
  nHp = nH - nH0
  
  y = Y / 4.0 / (1.0 - Y) # Eq. 32 in Katz & Weinberg - 1996 (Y = 0.24, i.e. He mass fraction. X = 1.0 - Y)
  nHepp = y * nH - nHe0 - nHep
  
  ne = nHp + nHep + 2.0 * nHepp
  
  dnH0_dt = k2(T) * nHp * ne - k1(T) * nH0 * ne
  
  dnHe0_dt = k5(T) * nHep * ne + k7(T) * nHep * ne - k3(T) * nHe0 * ne
  
  dnHep_dt = k6(T) * nHepp * ne + k3(T) * nHe0 * ne - k4(T) * nHep * ne - k5(T) * nHep * ne - k7(T) * nHep * ne
  
  Lamb = Lambda(T, nH, nH0, nHe0, nHep)
  
  ntot = n_tot(nH, nH0, nHe0, nHep)
  
  #------- Grassi et al - 2014 - eq. 8. - I excluded the nH2 (molecular) part. Add it if you consider molecules!!!!
  #nHe = nHe0 + nHep + nHepp
  #gamma = (5.*nH + 5.*nHe + 5.*ne) / (3.*nH + 3.*nHe + 3.*ne)
  #print(gamma)
  #-------

  dT_dt = -1.0 * (gamma - 1.0) / kB / ntot * Lamb
  
  return [dnH0_dt, dnHe0_dt, dnHep_dt, dT_dt]



gamma = 5./3.
kB = 1.3807e-16
X = 0.76
Y = 1.0 - X

nH = 1000.0

y0 = [1e-3, 1.6e-6, 6e-3, 1e6]

t_span = (1*3.16e7, 20000*3.16e7)

solution = solve_ivp(func, t_span, y0, method='LSODA', dense_output=True)

t = np.linspace(t_span[0], t_span[1], 10000) # This 10000 is not years, it is the number of points in linspace !!!!
y = solution.sol(t)

print(y.shape)

t_yrs = t / 3.16e7

nH0 = y[0, :]
nHe0 = y[1, :]
nHep = y[2, :]

yy = Y / 4.0 / (1.0 - Y) # Eq. 32 in Katz & Weinberg - 1996 (Y = 0.24, i.e. He mass fraction. X = 1.0 - Y)
nHepp = (yy * nH - nHe0 - nHep)

T = y[3, :]

nHp = (nH - nH0)

print('tot He = ', nHe0 + nHep + nHepp)

#----- Preparing cooling rate for plotting -----
res = []
for Tx, nH0x, nHe0x, nHepx in zip(T, nH0, nHe0, nHep):

  lmb = Lambda(Tx, nH, nH0x, nHe0x, nHepx)
  
  res.append([Tx, lmb])

res = np.array(res)

Tx = res[:, 0]
lmb = res[:, 1]


#------ Result from "test_primordial_hdf5_v2.py" code -----
with open('chimesRes.pkl', 'rb') as f:
  df = pickle.load(f)
# dictx = {'t_Arr_in_yrs': t_Arr_in_yrs, 'TEvol': TEvol, 'nHe0': nHe0, 'nHep': nHep, 'nHepp': nHepp}
t_Arr_in_yrsx = df['t_Arr_in_yrs']
TEvolx = df['TEvol']
nHe0x = df['nHe0']
nHepx = df['nHep']
nHeppx = df['nHepp']
#----------------------------------------------------------


plt.figure(figsize = (16, 14))

plt.subplot(2, 3, 1)
plt.scatter(t_yrs, np.log10(T), s = 5, color = 'k', label = 'my own code')
plt.scatter(t_Arr_in_yrsx, np.log10(TEvolx), s = 2, color = 'orange', label = 'chimes result', linestyle = '--')
plt.xlim(0, 20000)
plt.ylim(3, 8)
plt.legend()

plt.subplot(2, 3, 2)

plt.scatter(t_yrs, nHe0, s = 5, color = 'orange', label = 'nHe0')
plt.scatter(t_yrs, nHep, s = 5, color = 'lime', label = 'nHep')
plt.scatter(t_yrs, nHepp, s = 5, color = 'yellow', label = 'nHepp')
plt.yscale('log')
plt.title('solve_ivp')
plt.legend()


plt.subplot(2, 3, 3)
nHeTot = nHe0 + nHep + nHepp
plt.plot(T, nH0/nH, label = 'nH0', linestyle = ':')
plt.plot(T, nHp/nH, label = 'nHp', linestyle = ':')
plt.plot(T, nHe0/nHeTot, label = 'nHe0')
plt.plot(T, nHep/nHeTot, label = 'nHep')
plt.plot(T, nHepp/nHeTot,label = 'nHepp')
plt.yscale('log')
plt.xscale('log')
plt.ylim(2e-3, 1.2)
plt.xlim(1e4, 5e6)
plt.legend()

plt.subplot(2, 3, 4)
plt.scatter(np.log10(Tx), np.log10(lmb/nH/nH), s = 5, color = 'k')
plt.xlim(3.5, 8.25)
plt.ylim(-25, -21.5)

plt.tight_layout()

plt.savefig('myOnlyH.png')

plt.show()



