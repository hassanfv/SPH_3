
# In this version (i.e. mprog1.1.py) I used technique in Draine book section 27.3.1 for recombination cooling rates (Thanks to T. Grassi).
# This version (i.e. prog7.0.py) uses kB * T * k2(T) as the cooling due to recombination of Hp. We can convert recomb. rate to cooling rate by multiplying kB*T!
# This version (i.e. prog6.0.py) also includs cooling due to He atom and ions.
# This version also includes Free-Free cooling (So we have cooling due to Hydrogen atoms and free-free emission).
# This version uses solve_ivp.

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.integrate import odeint
import pickle


kB = 1.3807e-16 # erg/K
kB_in_eV_K = 8.61733326e-05 # eV/K


# Ref: Voronov - 1997 (parameters used for collisional ionization)
#                  dE    P      A        X      K
C_colIparams = [
                 [11.3, 0.0, 0.685e-7, 0.193, 0.25], # CI
                 [24.4, 1.0, 0.186e-7, 0.286, 0.24], # CII
                 [47.9, 1.0, 0.635e-8, 0.427, 0.21], # CIII
                 [64.5, 1.0, 0.150e-8, 0.416, 0.13], # CIV
                 [392.1,1.0, 0.299e-9, 0.666, 0.02], # CV
                 [490.0,1.0, 0.123e-9, 0.620, 0.16]  # CVI
               ]
          


#-----------------------------------------------------------

C_IP = np.array([11.3, 24.4, 47.9, 64.5, 392.1, 490.0]) * 1.60218e-12 # erg

# Ref: Shull & Steenberg - 1982. - Acol and Tcol are used for ionization BUT we use tables from Voronov - 1997 for collisional ionization rates!
# For Recombination       Acol      Tcol    Arad      Xrad     Adi      Bdi      T0      T1
C_RRparams = np.array([
                        [1.44e-10, 1.31e5, 4.70e-13, 6.24e-1, 2.54e-3, 4.42e-2, 1.57e5, 3.74e5], # CI
                        [4.20e-11, 2.83e5, 2.30e-12, 6.45e-1, 6.15e-3, 5.88e-2, 1.41e5, 1.41e5], # CII
                        [1.92e-11, 5.56e5, 3.20e-12, 7.70e-1, 1.62e-3, 3.43e-1, 8.19e4, 1.59e5], # CIII
                        [5.32e-12, 7.48e5, 7.50e-12, 8.17e-1, 4.78e-2, 3.62e-1, 3.44e6, 5.87e5], # CIV
                        [2.87e-13, 4.55e6, 1.70e-11, 7.21e-1, 3.22e-2, 3.15e-1, 4.06e6, 8.31e5], # CV
                        [9.16e-14, 5.68e6, 3.30e-11, 7.26e-1, 0.00e-0, 0.00e-0, 0.00e0, 0.00e0]  # CVI
                      ])


# Ref: Vonorov - 1997.
#----- phi_col (in cm^3.s^-1) (T in K, dE in eV) - It should be multiplied by ionization energy to become cooling rate!
# Z + e ---> Zp + 2e (collisional ionization)
def phi_col(T, i, j, params):

  dE, P, A, X, K = params[i]

  T_eV = kB_in_eV_K * T
  U = dE / T_eV

  return A * ((1. + P * U**0.5) / (X + U)) * U**K * np.exp(-U)


# Ref: Shull & Steenberg - 1982.
#----- RR_rate (Radiative Recombination rate)
def RR_rate(T, i, j, params):

  P = params[j]

  Arad = P[2]
  Xrad = P[3]
  
  Adi = P[4]
  Bdi = P[5]
  
  T0 = P[6]
  T1 = P[7]
  
  a_r = Arad / (T / 1e4)**Xrad
  a_d = Adi * T**(-3./2.) * np.exp(-T0/T) * (1. + Bdi * np.exp(-T1/T))

  return a_r + a_d


#----- getCarbonRates
def getCarbonRates(T, C_colIparams, C_RRparams):

  C_0_1 = phi_col(T, 0, 1, C_colIparams) # ionization of neutral Carbon CI to CII via collision:: CI + e ---> CII + 2e
  C_1_2 = phi_col(T, 1, 2, C_colIparams) # CII + e ---> CIII + 2e
  C_2_3 = phi_col(T, 2, 3, C_colIparams) # CIII + e ---> CIV + 2e
  C_3_4 = phi_col(T, 3, 4, C_colIparams) # CVI + e ---> CV + 2e
  C_4_5 = phi_col(T, 4, 5, C_colIparams) # CV + e ---> CVI + 2e
  C_5_6 = phi_col(T, 5, 6, C_colIparams) # CVI + e ---> CVII + 2e

  C_1_0 = RR_rate(T, 1, 0, C_RRparams) # recombination of singly ionized Carbon (CII) to neutral Carbon (CI):: CII + e ---> CI + photon
  C_2_1 = RR_rate(T, 2, 1, C_RRparams) # CIII + e ---> CII + photon
  C_3_2 = RR_rate(T, 3, 2, C_RRparams) # CIV + e ---> CIII + photon
  C_4_3 = RR_rate(T, 4, 3, C_RRparams) # CV + e ---> CIV + photon
  C_5_4 = RR_rate(T, 5, 4, C_RRparams) # CVI + e ---> CV + photon
  C_6_5 = RR_rate(T, 6, 5, C_RRparams) # CVII + e ---> CVI + photon

  return [
           C_0_1, C_1_2, C_2_3, C_3_4, C_4_5, C_5_6,
           C_1_0, C_2_1, C_3_2, C_4_3, C_5_4, C_6_5
         ]



#----- gamma_numeric
def gamma_numeric(T, i, j, params): # Only scalar T
  dT = T/10
  T1 = T - dT
  T2 = T + dT
  return (np.log(RR_rate(T2, i, j, params)) - np.log(RR_rate(T1, i, j, params))) / (np.log(T2) - np.log(T1)) - 0.5


#----- getKE_numeric
def getKE_numeric(T, i, j, params):
  gam = gamma_numeric(T, i, j, params)
  return (2.0 + gam) * kB * T


#----- getCarbonCoolingRates
def getCarbonCoolingRates(T, C_colIparams, C_RRparams):

  G_0_1 = phi_col(T, 0, 1, C_colIparams) * C_IP[0] # CI + e ---> CII + 2e
  G_1_2 = phi_col(T, 1, 2, C_colIparams) * C_IP[1] # CII + e ---> CIII + 2e
  G_2_3 = phi_col(T, 2, 3, C_colIparams) * C_IP[2] # CIII + e ---> CIV + 2e
  G_3_4 = phi_col(T, 3, 4, C_colIparams) * C_IP[3] # CVI + e ---> CV + 2e
  G_4_5 = phi_col(T, 4, 5, C_colIparams) * C_IP[4] # CV + e ---> CVI + 2e
  G_5_6 = phi_col(T, 5, 6, C_colIparams) * C_IP[5] # CVI + e ---> CVII + 2e


  G_1_0 = RR_rate(T, 1, 0, C_RRparams) * getKE_numeric(T, 1, 0, C_RRparams) # CII + e ---> CI + photon
  G_2_1 = RR_rate(T, 2, 1, C_RRparams) * getKE_numeric(T, 2, 1, C_RRparams) # CIII + e ---> CII + photon
  G_3_2 = RR_rate(T, 3, 2, C_RRparams) * getKE_numeric(T, 3, 2, C_RRparams) # CIV + e ---> CIII + photon
  G_4_3 = RR_rate(T, 4, 3, C_RRparams) * getKE_numeric(T, 4, 3, C_RRparams) # CV + e ---> CIV + photon
  G_5_4 = RR_rate(T, 5, 4, C_RRparams) * getKE_numeric(T, 5, 4, C_RRparams) # CVI + e ---> CV + photon
  G_6_5 = RR_rate(T, 6, 5, C_RRparams) * getKE_numeric(T, 6, 5, C_RRparams) # CVII + e ---> CVI + photon

  return [
           G_0_1, G_1_2, G_2_3, G_3_4, G_4_5, G_5_6,
           G_1_0, G_2_1, G_3_2, G_4_3, G_5_4, G_6_5
         ]




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
def Lambda(T, nH, nH0, nHe0, nHep, nC, nC0, nC1, nC2, nC3, nC4, nC5):

  nHp = nH - nH0
  y = Y / 4.0 / (1.0 - Y) # Eq. 32 in Katz & Weinberg - 1996 (Y = 0.24, i.e. He mass fraction. X = 1.0 - Y)
  nHepp = y * nH - nHe0 - nHep
  
  nC6 = nC - (nC0 + nC1 + nC2 + nC3 + nC4 + nC5)
  
  ne = nHp + (nHep + 2.0 * nHepp) + (nC1 + 2.0 * nC2 + 3.0 * nC3 + 4.0 * nC4 + 5.0 * nC5 + 6.0 * nC6)
  
  C_0_1, C_1_2, C_2_3, C_3_4, C_4_5, C_5_6, C_1_0, C_2_1, C_3_2, C_4_3, C_5_4, C_6_5  = getCarbonRates(T, C_colIparams, C_RRparams)
  G_0_1, G_1_2, G_2_3, G_3_4, G_4_5, G_5_6, G_1_0, G_2_1, G_3_2, G_4_3, G_5_4, G_6_5  = getCarbonCoolingRates(T, C_colIparams, C_RRparams)
  
  ff_coeff = nC1 + 2**2*nC2 + 3**2*nC3 + 4**2*nC4 + 5**2*nC5 + 6**2*nC6
  
  Lamb = (
           g1(T) * ne * nH0 # collisional ionization of hydrogen ::::::: (H0 + e ---> Hp + 2e).
         + g2(T)  * ne * nHp # photo-recombination of Hp with electron :: (Hp + e ---> H0 + γ).
         + g3(T)  * ne * nH0 # collisional excitaion of H0.
         + g4(T) * ne * (nHp + nHep + 4.0 * nHepp + ff_coeff) # free-free emission !!!!!!!!! How to incorporate metals contribution??????????
         + g5(T) * nHep * ne # collisional excitation of Hep.
         + g6(T) * nHe0 * ne # He0 collisional ionization
         + g7(T) * nHep * ne # Hep collisional ionization
         + g8(T) * nHep * ne # Hep recombination to He0
         + g9(T) * nHepp * ne# Hepp recombination to Hep
         + g10(T) * nHep * ne# Hep di-electric recombination to He0
         + G_0_1 * ne * nC0 + G_1_0 * ne * nC1
         + G_1_2 * ne * nC1 + G_2_1 * ne * nC2
         + G_2_3 * ne * nC2 + G_3_2 * ne * nC3
         + G_3_4 * ne * nC3 + G_4_3 * ne * nC4
         + G_4_5 * ne * nC4 + G_5_4 * ne * nC5
         + G_5_6 * ne * nC5 + G_6_5 * ne * nC6
         )
  return Lamb



#----- n_tot
def n_tot(nH, nH0, nHe0, nHep, nC, nC0, nC1, nC2, nC3, nC4, nC5):
  
  nHp = nH - nH0

  y = Y / 4.0 / (1.0 - Y) # Eq. 32 in Katz & Weinberg - 1996 (Y = 0.24, i.e. He mass fraction. X = 1.0 - Y)
  nHepp = y * nH - nHe0 - nHep
  
  nC6 = nC - (nC0 + nC1 + nC2 + nC3 + nC4 + nC5)
  
  ne = nHp + (nHep + 2.0 * nHepp) + (nC1 + 2.0 * nC2 + 3.0 * nC3 + 4.0 * nC4 + 5.0 * nC5 + 6.0 * nC6)
  
  ntot = ne + (nH0 + nHp) + (nHe0 + nHep + nHepp) + (nC0 + nC1 + nC2 + nC3 + nC4 + nC5 + nC6)
  
  return ntot



#-----------------------------------------
#------ Solution of the ODE Section ------
#-----------------------------------------
def func(t, y):

  nH0, nHe0, nHep, nC0, nC1, nC2, nC3, nC4, nC5, T = y
  
  nHp = nH - nH0
  
  nC6 = nC - (nC0 + nC1 + nC2 + nC3 + nC4 + nC5)
  
  y = Y / 4.0 / (1.0 - Y) # Eq. 32 in Katz & Weinberg - 1996 (Y = 0.24, i.e. He mass fraction. X = 1.0 - Y)
  nHepp = y * nH - nHe0 - nHep
  
  ne = nHp + (nHep + 2.0 * nHepp) + (nC1 + 2.0 * nC2 + 3.0 * nC3 + 4.0 * nC4 + 5.0 * nC5 + 6.0 * nC6)
  
  C_0_1, C_1_2, C_2_3, C_3_4, C_4_5, C_5_6, C_1_0, C_2_1, C_3_2, C_4_3, C_5_4, C_6_5  = getCarbonRates(T, C_colIparams, C_RRparams)
  
  dnH0_dt = k2(T) * nHp * ne - k1(T) * nH0 * ne
  dnHe0_dt = k5(T) * nHep * ne + k7(T) * nHep * ne - k3(T) * nHe0 * ne
  dnHep_dt = k6(T) * nHepp * ne + k3(T) * nHe0 * ne - k4(T) * nHep * ne - k5(T) * nHep * ne - k7(T) * nHep * ne
  
  dnC0_dt = C_1_0 * ne * nC1 - C_0_1 * ne * nC0
  dnC1_dt = C_0_1 * ne * nC0 + C_2_1 * ne * nC2 - C_1_0 * ne * nC1 - C_1_2 * ne * nC1
  dnC2_dt = C_1_2 * ne * nC1 + C_3_2 * ne * nC3 - C_2_1 * ne * nC2 - C_2_3 * ne * nC2
  dnC3_dt = C_2_3 * ne * nC2 + C_4_3 * ne * nC4 - C_3_2 * ne * nC3 - C_3_4 * ne * nC3
  dnC4_dt = C_3_4 * ne * nC3 + C_5_4 * ne * nC5 - C_4_3 * ne * nC4 - C_4_5 * ne * nC4
  dnC5_dt = C_4_5 * ne * nC4 + C_6_5 * ne * nC6 - C_5_4 * ne * nC5 - C_5_6 * ne * nC5
  
  Lamb = Lambda(T, nH, nH0, nHe0, nHep, nC, nC0, nC1, nC2, nC3, nC4, nC5)
  
  ntot = n_tot(nH, nH0, nHe0, nHep, nC, nC0, nC1, nC2, nC3, nC4, nC5)

  dT_dt = -1.0 * (gamma - 1.0) / kB / ntot * Lamb
  
  return [dnH0_dt, dnHe0_dt, dnHep_dt, dnC0_dt, dnC1_dt, dnC2_dt, dnC3_dt, dnC4_dt, dnC5_dt, dT_dt]








gamma = 5./3.
kB = 1.3807e-16
X = 0.76
Y = 1.0 - X

nH = 1000.0

C_solar = 10**(-3.57)
nC = C_solar * nH

print('nC (cm^-3) = ', nC)

print()
print('H, C before = ', nH, nC)

#      nH0   nHe0   nHep   nC0  nC1    nC2   nC3   nC4  nC5    T
y0 = [1e-4, 1.3e-8, 4e-4, 1e-5, 1e-5, 1e-5, 1e-2, 1e-2, 1e-2, 1e6]

t_span = (1*3.16e7, 20000*3.16e7)

solution = solve_ivp(func, t_span, y0, method='LSODA', dense_output=True)

t = np.linspace(t_span[0], t_span[1], 10000)
y = solution.sol(t)

print(y.shape)

t_yrs = t / 3.16e7

nH0 = y[0, :]
nHe0 = y[1, :]
nHep = y[2, :]
nC0 = y[3, :]
nC1 = y[4, :]
nC2 = y[5, :]
nC3 = y[6, :]
nC4 = y[7, :]
nC5 = y[8, :]
nC6 = nC - (nC0 + nC1 + nC2 + nC3 + nC4 + nC5)
T = y[9, :]

yy = Y / 4.0 / (1.0 - Y) # Eq. 32 in Katz & Weinberg - 1996 (Y = 0.24, i.e. He mass fraction. X = 1.0 - Y)
nHepp = (yy * nH - nHe0 - nHep)

nHp = (nH - nH0)


#----- Preparing cooling rate for plotting -----

res = []
for Tx, nH0x, nHe0x, nHepx, nC0x, nC1x, nC2x, nC3x, nC4x, nC5x in zip(T, nH0, nHe0, nHep, nC0, nC1, nC2, nC3, nC4, nC5):

  lmb = Lambda(Tx, nH, nH0x, nHe0x, nHepx, nC, nC0x, nC1x, nC2x, nC3x, nC4x, nC5x)
  
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
nH0x = df['nH0']
nHpx = df['nHp']
nHx = nH0x + nHpx
nHe0x = df['nHe0']
nHepx = df['nHep']
nHeppx = df['nHepp']
nHeTotx = nHe0x + nHepx + nHeppx

nC0x = df['nC0']
nC1x = df['nC1']
nC2x = df['nC2']
nC3x = df['nC3']
nC4x = df['nC4']
nC5x = df['nC5']
nC6x = df['nC6']
nCx = nC0x + nC1x + nC2x + nC3x + nC4x + nC5x + nC6x
#----------------------------------------------------------


plt.figure(figsize = (16, 8))

plt.subplot(2, 3, 1)
plt.scatter(t_yrs, np.log10(T), s = 5, color = 'k', label = 'my own code')
plt.scatter(t_Arr_in_yrsx, np.log10(TEvolx), s = 2, color = 'orange', label = 'chimes result', linestyle = '--')
plt.xlim(0, 3000)
plt.ylim(3, 8)
plt.legend()

plt.subplot(2, 3, 2)

plt.plot(t_yrs, nHe0, color = 'r', label = 'nHe0')
plt.plot(t_yrs, nHep, color = 'g', label = 'nHep')
plt.plot(t_yrs, nHepp, color = 'b', label = 'nHepp')

plt.plot(t_Arr_in_yrsx, nHe0x, color = 'r', label = 'nHe0 - chimes', linestyle = ':')
plt.plot(t_Arr_in_yrsx, nHepx, color = 'g', label = 'nHep - chimes', linestyle = ':')
plt.plot(t_Arr_in_yrsx, nHeppx, color = 'b', label = 'nHepp - chimes', linestyle = ':')

plt.xlim(0, 5000)

plt.yscale('log')
plt.title('solve_ivp')
plt.legend()


plt.subplot(2, 3, 3)
nHeTot = nHe0 + nHep + nHepp
plt.plot(T, nH0/nH, label = 'nH0', color = 'r')
plt.plot(T, nHp/nH, label = 'nHp', color = 'g')
plt.plot(T, nHe0/nHeTot, label = 'nHe0', color = 'b')
plt.plot(T, nHep/nHeTot, label = 'nHep', color = 'orange')
plt.plot(T, nHepp/nHeTot,label = 'nHepp', color = 'purple')

plt.plot(TEvolx, nH0x/nHx, label = 'nH0 - chimes', color = 'r', linestyle = ':')
plt.plot(TEvolx, nHpx/nHx, label = 'nHp - chimes', color = 'g', linestyle = ':')
plt.plot(TEvolx, nHe0x/nHeTotx, label = 'nHe0 - chimes', color = 'b', linestyle = ':')
plt.plot(TEvolx, nHepx/nHeTotx, label = 'nHep - chimes', color = 'orange', linestyle = ':')
plt.plot(TEvolx, nHeppx/nHeTotx,label = 'nHepp - chimes', color = 'purple', linestyle = ':')

plt.yscale('log')
plt.xscale('log')
plt.ylim(2e-3, 1.2)
plt.xlim(1e4, 1e6)
plt.legend()

plt.subplot(2, 3, 4)
plt.scatter(np.log10(Tx), np.log10(lmb/nH/nH), s = 5, color = 'k')
plt.xlim(3.5, 8.25)
plt.ylim(-25, -21.5)

ne = nHp + (nHep + 2.0 * nHepp) + (nC1 + 2.0 * nC2 + 3.0 * nC3 + 4.0 * nC4 + 5.0 * nC5 + 6.0 * nC6)
#nCC = (nC0 + nC1 + nC2 + nC3 + nC4 + nC5 + nC6)
print('nH after = ', nH)
print('ne after (one ne for each T) = ', ne)
#print('nC after (one nC for each T) = ', nCC)

plt.subplot(2, 3, 5)
plt.plot(T, nC0/nC, label = 'nC0', color = 'r')
plt.plot(T, nC1/nC, label = 'nC1', color = 'g')
plt.plot(T, nC2/nC, label = 'nC2', color = 'b')
plt.plot(T, nC3/nC, label = 'nC3', color = 'orange')
plt.plot(T, nC4/nC, label = 'nC4', color = 'purple')
plt.plot(T, nC5/nC, label = 'nC5', color = 'lime')
plt.plot(T, nC6/nC, label = 'nC6', color = 'pink')

plt.plot(TEvolx, nC0x/nCx, color = 'r', linestyle = ':')
plt.plot(TEvolx, nC1x/nCx, color = 'g', linestyle = ':')
plt.plot(TEvolx, nC2x/nCx, color = 'b', linestyle = ':')
plt.plot(TEvolx, nC3x/nCx, color = 'orange', linestyle = ':')
plt.plot(TEvolx, nC4x/nCx, color = 'purple', linestyle = ':')
plt.plot(TEvolx, nC5x/nCx, color = 'lime', linestyle = ':')
plt.plot(TEvolx, nC6x/nCx, color = 'pink', linestyle = ':')

plt.yscale('log')
plt.xscale('log')
plt.ylim(2e-3, 2.2)
plt.xlim(1e4, 3e6)
plt.legend()

plt.tight_layout()

plt.savefig('myH_He_C.png')

plt.show()



