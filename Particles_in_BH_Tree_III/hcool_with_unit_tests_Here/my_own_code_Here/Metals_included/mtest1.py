
# Ref: Grassi et al. - 2011 - ROBO code!
# Ref for recombination rates: Table 3 in Grassi et al. - 2011 or Verner & Ferland - 1996.

import numpy as np
import matplotlib.pyplot as plt


# Cooling rate due to H collisional ionization: (H0 + e ---> Hp + 2e)  ::: in erg.s^-1.cm^3  ::: will be multiplied by nH0 and ne later in the code.
def g1(T):

  g1_val = 1.27e-21 * T**0.5 * np.exp(-157809.1/T) / (1.0 + (T/1e5)**0.5)
  
  return g1_val


# Cooling rate due to Hp photo-recombination: (Hp + e ---> H0 + γ)  ::: in erg.s^-1.cm^3  ::: will be multiplied by nHp and ne later in the code.
def g2(T):

  g2_val = 8.70e-27 * T**0.5 / (T/1e3)**0.2 / (1.0 + (T/1e6)**0.7)
  #g2_val = kB * T * k2(T) # See LaMothe and Ferland - 2001, page 1, the lines below equation 1!

  return g2_val




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



CVI_IP = 489.99334 * 1.60218e-12 # erg. IP of ionizaing CV to CVI (i.e. fully ionized Carbon).

C_rec = [6.556e-10, 65.23, 2.446e7, 0.7567]
C_col = [6.85e-8, 0, 11.3, 0.193, 0.25]

Tgrid = np.logspace(4, 8, 100)

res = []
for T in Tgrid:
  res.append([T, kB * T * phi_rec(T, C_rec), CVI_IP * phi_col(T, C_col), g2(T), g1(T)])

res = np.array(res)

T = res[:, 0]
C_rec_rate = res[:, 1]
C_col_rate = res[:, 2]

H_rec_rate = res[:, 3]
H_col_rate = res[:, 4]

plt.scatter(np.log10(T), np.log10(C_rec_rate), s = 5, color = 'k', label = 'C_rec_rate')
plt.scatter(np.log10(T), np.log10(C_col_rate), s = 5, color = 'b', label = 'C_col_rate')

plt.scatter(np.log10(T), np.log10(H_rec_rate), s = 5, color = 'cyan', label = 'H_rec_rate')
plt.scatter(np.log10(T), np.log10(H_col_rate), s = 5, color = 'r', label = 'H_col_rate')
plt.legend()

plt.show()




