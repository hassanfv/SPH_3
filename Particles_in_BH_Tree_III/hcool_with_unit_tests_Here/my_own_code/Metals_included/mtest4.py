
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
  C_3_4 = C_IP[3] * phi_col(T, 3, 4, C_colIparams) # CVI + e ---> CV + 2e
  C_4_5 = C_IP[4] * phi_col(T, 4, 5, C_colIparams) # CV + e ---> CVI + 2e
  C_5_6 = C_IP[5] * phi_col(T, 5, 6, C_colIparams) # CVI + e ---> CVII + 2e

  C_1_0 = kB * T * RR_rate(T, 1, 0, C_RRparams) # recombination of singly ionized Carbon (CII) to neutral Carbon (CI):: CII + e ---> CI + photon
  C_2_1 = kB * T * RR_rate(T, 2, 1, C_RRparams) # CIII + e ---> CII + photon
  C_3_2 = kB * T * RR_rate(T, 3, 2, C_RRparams) # CIV + e ---> CIII + photon
  C_4_3 = kB * T * RR_rate(T, 4, 3, C_RRparams) # CV + e ---> CIV + photon
  C_5_4 = kB * T * RR_rate(T, 5, 4, C_RRparams) # CVI + e ---> CV + photon
  C_6_5 = kB * T * RR_rate(T, 6, 5, C_RRparams) # CVII + e ---> CVI + photon

  return [
           C_0_1, C_1_2, C_2_3, C_3_4, C_4_5, C_5_6,
           C_1_0, C_2_1, C_3_2, C_4_3, C_5_4, C_6_5
         ]







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




Tgrid = np.logspace(3, 8, 100)

res = []
for T in Tgrid:

  C_0_1, C_1_2, C_2_3, C_3_4, C_4_5, C_5_6, C_1_0, C_2_1, C_3_2, C_4_3, C_5_4, C_6_5  = getCarbonRates(T, C_colIparams, C_RRparams)

  res.append([C_0_1, C_1_2, C_2_3, C_3_4, C_4_5, C_5_6, C_1_0, C_2_1, C_3_2, C_4_3, C_5_4, C_6_5])


res = np.array(res)

C_0_1 = res[:, 0] # Collisional Ionization rate for: CI + e ---> CII + 2e
C_1_2 = res[:, 1] # Collisional Ionization rate for: CII + e ---> CIII + 2e
C_2_3 = res[:, 2] # Collisional Ionization rate for: CIII + e ---> CIV + 2e
C_3_4 = res[:, 3] # Collisional Ionization rate for: CIV + e ---> CV + 2e
C_4_5 = res[:, 4] # Collisional Ionization rate for: CV + e ---> CVI + 2e
C_5_6 = res[:, 5] # Collisional Ionization rate for: CVI + e ---> CVII + 2e

C_1_0 = res[:, 6]  # Recombination rate for : CII + e ---> CI + photon
C_2_1 = res[:, 7]  # Recombination rate for : CIII + e ---> CII + photon
C_3_2 = res[:, 8]  # Recombination rate for : CIV + e ---> CIII + photon
C_4_3 = res[:, 9]  # Recombination rate for : CV + e ---> CIV + photon
C_5_4 = res[:, 10] # Recombination rate for : CVI + e ---> CV + photon
C_6_5 = res[:, 11] # Recombination rate for : CVII + e ---> CVI + photon


plt.plot(np.log10(Tgrid), np.log10(C_0_1), color = 'r', label = 'C_0_1', linestyle = '--')
plt.plot(np.log10(Tgrid), np.log10(C_1_2), color = 'g', label = 'C_1_2', linestyle = '--')
plt.plot(np.log10(Tgrid), np.log10(C_2_3), color = 'b', label = 'C_2_3', linestyle = '--')
plt.plot(np.log10(Tgrid), np.log10(C_3_4), color = 'cyan', label = 'C_3_4', linestyle = '--')
plt.plot(np.log10(Tgrid), np.log10(C_4_5), color = 'lime', label = 'C_4_5', linestyle = '--')
plt.plot(np.log10(Tgrid), np.log10(C_5_6), color = 'purple', label = 'C_5_6', linestyle = '--')

plt.plot(np.log10(Tgrid), np.log10(C_1_0), color = 'r', label = 'C_1_0')
plt.plot(np.log10(Tgrid), np.log10(C_2_1), color = 'g', label = 'C_2_1')
plt.plot(np.log10(Tgrid), np.log10(C_3_2), color = 'b', label = 'C_3_2')
plt.plot(np.log10(Tgrid), np.log10(C_4_3), color = 'cyan', label = 'C_4_3')
plt.plot(np.log10(Tgrid), np.log10(C_5_4), color = 'lime', label = 'C_5_4')
plt.plot(np.log10(Tgrid), np.log10(C_6_5), color = 'purple', label = 'C_6_5')

plt.ylim(-25, -17)

plt.legend()

plt.savefig('fig.png')

plt.show()




