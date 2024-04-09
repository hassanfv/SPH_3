
# Ref: Grassi et al. - 2011 - ROBO code!
# Ref for recombination rates: Table 3 in Grassi et al. - 2011 or Verner & Ferland - 1996.

import numpy as np
import matplotlib.pyplot as plt


kB = 1.3807e-16 # erg/K
kB_in_eV_K = 8.61733326e-05 # eV/K

#----- ColI_rate (Collisional Ionization rate) NOTE: I just have this here for comparison with Voronov - 1997 rates (as unit test!)!
def ColI_rate(T, P):
  
  Acol = P[0]
  Tcol = P[1]

  ai = 0.1
  
  return Acol * T**0.5 / (1.0 + ai * T / Tcol) * np.exp(-Tcol/T)

#----- phi_rec (in cm^3.s^-1) It should be multiplied by kB*T to become cooling rate! Used just for testing!
def phi_rec(T, Zrec):

  A, tou0, tou1, b = Zrec

  T0 = (T/tou0)**0.5
  T1 = (T/tou1)**0.5

  return A / (T0 * (1. + T0)**(1. - b) * (1. + T1)**(1. + b))
#---------------------------------------------------------------------------------------


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
  
  a_r = Arad * (T / 1e4)**-Xrad
  a_d = Adi * T**-(3./2.) * np.exp(-T0/T) * (1. + Bdi * np.exp(-T1/T))

  return a_r + a_d


#----- getCarbonRates
def getCarbonRates(T, C_colIparams, C_RRparams):

  C_0_1 = C_IP[0] * phi_col(T, 0, 1, C_colIparams) # ionization of neutral Carbon CI to CII via collision:: CI + e ---> CII + 2e
  C_1_2 = C_IP[1] * phi_col(T, 1, 2, C_colIparams) # CII + e ---> CIII + 2e
  C_2_3 = C_IP[2] * phi_col(T, 2, 3, C_colIparams) # CIII + e ---> CIV + 2e

  C_1_0 = kB * T * RR_rate(T, 1, 0, C_RRparams) # recombination of singly ionized Carbon (CII) to neutral Carbon (CI):: CII + e ---> CI + photon
  C_2_1 = kB * T * RR_rate(T, 2, 1, C_RRparams) # CIII + e ---> CII + photon
  C_3_2 = kB * T * RR_rate(T, 3, 2, C_RRparams) # CIV + e ---> CIII + photon

  return [
           C_0_1, C_1_2, C_2_3,
           C_1_0, C_2_1, C_3_2
         ]







# Ref: Voronov - 1997
C_colIparams = [
          [11.3, 0.0, 0.685e-7, 0.193, 0.25],
          [24.4, 1.0, 0.186e-7, 0.286, 0.24],
          [47.9, 1.0, 0.635e-8, 0.427, 0.21]
        ]
          


#-----------------------------------------------------------

C_IP = np.array([11.3, 24.4, 47.9]) * 1.60218e-12 # erg

# Ref: Shull & Steenberg - 1982. - Acol and Tcol are used for ionization BUT we use tables from Voronov - 1997 for collisional ionization rates!
#                     Acol      Tcol    Arad      Xrad     Adi      Bdi      T0      T1
C_RRparams = np.array([
                     [1.44e-10, 1.31e5, 4.70e-13, 6.24e-1, 2.54e-3, 4.42e-2, 1.57e5, 3.74e5],
                     [4.20e-11, 2.83e5, 2.30e-12, 6.45e-1, 6.15e-3, 5.88e-2, 1.41e5, 1.41e5],
                     [1.92e-11, 5.56e5, 3.20e-12, 7.70e-1, 1.62e-3, 3.43e-1, 8.19e4, 1.59e5]
                   ])




Tgrid = np.logspace(3, 8, 100)

res = []
for T in Tgrid:

  C_0_1, C_1_2, C_2_3, C_1_0, C_2_1, C_3_2  = getCarbonRates(T, C_colIparams, C_RRparams)

  res.append([C_0_1, C_1_2, C_2_3, C_1_0, C_2_1, C_3_2])


res = np.array(res)

C_0_1 = res[:, 0] # Collisional Ionization rate for: CI + e ---> CII + 2e
C_1_2 = res[:, 1] # Collisional Ionization rate for: CII + e ---> CIII + 2e
C_2_3 = res[:, 2] # Collisional Ionization rate for: CIII + e ---> CIV + 2e

C_1_0 = res[:, 3] # Recombination rate for : CII + e ---> CI + photon
C_2_1 = res[:, 4] # Recombination rate for : CIII + e ---> CII + photon
C_3_2 = res[:, 5] # Recombination rate for : CIV + e ---> CIII + photon


plt.plot(np.log10(Tgrid), np.log10(C_0_1), color = 'red', label = 'C_0_1')
plt.plot(np.log10(Tgrid), np.log10(C_1_2), color = 'red', label = 'C_1_2', linestyle = '--')
plt.plot(np.log10(Tgrid), np.log10(C_2_3), color = 'red', label = 'C_2_3', linestyle = ':')

plt.plot(np.log10(Tgrid), np.log10(C_1_0), color = 'blue', label = 'C_1_0')
plt.plot(np.log10(Tgrid), np.log10(C_2_1), color = 'blue', label = 'C_2_1', linestyle = '--')
plt.plot(np.log10(Tgrid), np.log10(C_3_2), color = 'blue', label = 'C_3_2', linestyle = ':')

plt.ylim(-25, -17)

plt.legend()

plt.show()




