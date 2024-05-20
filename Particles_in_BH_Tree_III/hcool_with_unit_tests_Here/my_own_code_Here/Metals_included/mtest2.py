
# Ref: Grassi et al. - 2011 - ROBO code!
# Ref for recombination rates: Table 3 in Grassi et al. - 2011 or Verner & Ferland - 1996.

import numpy as np
import matplotlib.pyplot as plt


kB = 1.3807e-16 # erg/K
kB_in_eV_K = 8.61733326e-05 # eV/K

#----- phi_col (in cm^3.s^-1) (T in K, dE in eV) - It should be multiplied by ionization energy to become cooling rate!
# Z + e ---> Zp + 2e (collisional ionization)
def phi_col(T, Zcol):

  A, P, dE, X, K = Zcol

  T_eV = kB_in_eV_K * T
  U = dE / T_eV

  return A * ((1. + P * U**0.5) / (X + U)) * U**K * np.exp(-U)


#----- phi_rec (in cm^3.s^-1) It should be multiplied by kB*T to become cooling rate!
def phi_rec(T, Zrec):

  A, tou0, tou1, b = Zrec

  T0 = (T/tou0)**0.5
  T1 = (T/tou1)**0.5

  return A / (T0 * (1. + T0)**(1. - b) * (1. + T1)**(1. + b))



CI_IP = 11.26030 * 1.60218e-12 # erg. IP of ionizaing CV to CVI (i.e. fully ionized Carbon).

C_rec = [6.556e-10, 65.23, 2.446e7, 0.7567]
C_col = [6.85e-8, 0, 11.3, 0.193, 0.25]

Tgrid = np.logspace(3, 8, 100)

res = []
for T in Tgrid:
  res.append([T, kB * T * phi_rec(T, C_rec), CI_IP * phi_col(T, C_col)])

res = np.array(res)

T = res[:, 0]
C_rec_rate = res[:, 1]
C_col_rate = res[:, 2]


#-----------------------------------------------------------
Cparams = np.array([[1.44e-10, 1.31e5, 4.70e-13, 6.24e-1, 2.54e-3, 4.42e-2, 1.57e5, 3.74e5],
                    [4.20e-11, 2.83e5, 2.30e-12, 6.45e-1, 6.15e-3, 5.88e-2, 1.41e5, 1.41e5]
                    ])


#----- RR_rate (Radiative Recombination rate)
def RR_rate(T, P):

  Arad = P[2]
  Xrad = P[3]
  
  Adi = P[4]
  Bdi = P[5]
  
  T0 = P[6]
  T1 = P[7]
  
  a_r = Arad * (T / 1e4)**-Xrad
  a_d = Adi * T**-(3./2.) * np.exp(-T0/T) * (1. + Bdi * np.exp(-T1/T))

  return a_r + a_d


#----- ColI_rate (Collisional Ionization rate) NOTE:
def ColI_rate(T, P):
  
  Acol = P[0]
  Tcol = P[1]

  ai = 0.1
  
  return Acol * T**0.5 / (1.0 + ai * T / Tcol) * np.exp(-Tcol/T)



C1_IP = 11.26030 * 1.60218e-12 # erg ---> ionozation potential energy of CI (energy needed to convert CI to CII)
C2_IP = 24.38332 * 1.60218e-12 # erg ---> ionozation potential energy of CII (energy needed to convert CII to CIII)


res = []
for Tx in Tgrid:
  C1RR_rate = kB * Tx * RR_rate(Tx, Cparams[0])
  C2RR_rate = kB * Tx * RR_rate(Tx, Cparams[1])
  
  C1CI_rate = C1_IP * ColI_rate(Tx, Cparams[0])
  C2CI_rate = C2_IP * ColI_rate(Tx, Cparams[1])

  res.append([C1RR_rate, C2RR_rate, C1CI_rate, C2CI_rate])


res = np.array(res)

C1RR = res[:, 0] # Recombination rate for : CII + e ---> CI + photon
C2RR = res[:, 1] # Recombination rate for : CIII + e ---> CII + photon

C1ColI = res[:, 2] # Collisional Ionization rate for: CI + e ---> CII + 2e
C2ColI = res[:, 3] # Collisional Ionization rate for: CII + e ---> CIII + 2e

plt.scatter(np.log10(T), np.log10(C_rec_rate), s = 5, color = 'k', label = 'C_rec_rate')
plt.scatter(np.log10(T), np.log10(C_col_rate), s = 5, color = 'b', label = 'C_col_rate')

plt.scatter(np.log10(T), np.log10(C1RR), s = 5, color = 'orange', label = 'C1_RR New')
plt.scatter(np.log10(T), np.log10(C2RR), s = 5, color = 'red', label = 'C2_RR New')

plt.scatter(np.log10(T), np.log10(C1ColI), s = 5, color = 'yellow', label = 'C1_ColI New')
plt.scatter(np.log10(T), np.log10(C2ColI), s = 5, color = 'lime', label = 'C2_ColI New')

plt.ylim(-25, -17)

plt.legend()

plt.show()




