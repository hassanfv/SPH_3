
#Ref: https://www.youtube.com/watch?v=vNoFdtcPFdk

import numpy as np
import matplotlib.pyplot as plt
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
  nHepp = nHe - nHe0 - nHep
  
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

  nHepp = nHe - nHe0 - nHep
  
  ne = nHp + nHep + 2.0 * nHepp
  
  ntot = nH0 + nHp + nHe0 + nHep + nHepp + ne
  
  return ntot



#===== rk4SingleStep
def rk4SingleStep(func, dt, t0, y0):

  f1 = func(t0, y0)
  f2 = func(t0 + dt/2.0, y0 + (dt/2.0) * f1)
  f3 = func(t0 + dt/2.0, y0 + (dt/2.0) * f2)
  f4 = func(t0 + dt, y0 + dt * f3)
  
  yout = y0 + (dt / 6.0) * (f1 + 2.0 * f2 + 2.0 * f3 + f4)
  
  return yout



#-----------------------------------------
#------ Solution of the ODE Section ------
#-----------------------------------------
def func(t, y):

  nH0, nHe0, nHep, T = y
  
  nHp = nH - nH0

  nHepp = nHe - nHe0 - nHep
  
  ne = nHp + nHep + 2.0 * nHepp
  
  dnH0_dt = k2(T) * nHp * ne - k1(T) * nH0 * ne
  
  dnHe0_dt = k5(T) * nHep * ne + k7(T) * nHep * ne - k3(T) * nHe0 * ne
  
  dnHep_dt = k6(T) * nHepp * ne + k3(T) * nHe0 * ne - k4(T) * nHep * ne - k5(T) * nHep * ne - k7(T) * nHep * ne
  
  Lamb = Lambda(T, nH, nH0, nHe0, nHep)
  
  ntot = n_tot(nH, nH0, nHe0, nHep)

  dT_dt = -1.0 * (gamma - 1.0) / kB / ntot * Lamb
  
  return np.array([dnH0_dt, dnHe0_dt, dnHep_dt, dT_dt])




#----- p_interp
def g_interp(Tref, T0, gx):

  nt = -1
  
  N = len(Tref)
  
  for j in range(N):
    
    if (Tref[j] >= T0) & (nt == -1):
      nt = j
  
  if nt == -1:
    nt = N - 1 
  
  diff = Tref[nt] - Tref[nt - 1]
  fhi = (T0 - Tref[nt - 1]) / diff
  flow = 1.0 - fhi
  
  gnew = flow * np.log10(gx[nt - 1]) + fhi * np.log10(gx[nt])
  
  return 10**gnew
  




gamma = 5./3.
kB = 1.3807e-16
X = 0.76
Y = 1.0 - X

nH = 1000.0

Tref = np.logspace(4, 8, 100)

g1x = np.zeros(len(Tref))
g2x = np.zeros(len(Tref))
g3x = np.zeros(len(Tref))
g4x = np.zeros(len(Tref))
g5x = np.zeros(len(Tref))
g6x = np.zeros(len(Tref))
g7x = np.zeros(len(Tref))
g8x = np.zeros(len(Tref))
g9x = np.zeros(len(Tref))
g10x = np.zeros(len(Tref))

for i in range(len(Tref)):
  g1x[i] = g1(Tref[i])
  g2x[i] = g2(Tref[i])
  g3x[i] = g3(Tref[i])
  g4x[i] = g4(Tref[i])
  g5x[i] = g5(Tref[i])
  g6x[i] = g6(Tref[i])
  g7x[i] = g7(Tref[i])
  g8x[i] = g8(Tref[i])
  g9x[i] = g9(Tref[i])
  g10x[i] = g10(Tref[i])


T0 = 11213.0
g1new = g_interp(Tref, T0, g1x)

print(T0, g1(T0), g1new)

dictx = {'Tref': Tref, 'g1': g1x, 'g2': g2x, 'g3': g3x, 'g4': g4x, 'g5': g5x, 'g6': g6x, 'g7': g7x, 'g8': g8x, 'g9': g9x, 'g10': g10x}

with open('HeH_cooling_rates.pkl', 'wb') as f:
  pickle.dump(dictx, f)







