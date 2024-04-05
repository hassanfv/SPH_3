
# Ref: Grassi et al - 2011

import numpy as np
import matplotlib.pyplot as plt


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
  
  k1_val = np.exp(-32.71396786
                  + 13.5365560 * ln_T
                  - 5.73932875 * ln_T**2
                  + 1.56315498 * ln_T**3
                  - 0.287705600 * ln_T**4
                  + 3.48255977 * 10**-2 * ln_T**5
                  - 2.63197617 * 10**-3 * ln_T**6
                  + 1.11954395 * 10**-4 * ln_T**7
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
  
  g1_val = 13.5984 * 1.60218e-12 * k1
  
  return g1_val


# Cooling rate due to Hp photo-recombination: (Hp + e ---> H0 + γ)  ::: in erg.s^-1.cm^3  ::: will be multiplied by nHp and ne later in the code.
def g2(T):

  g2_val = 13.5984 * 1.60218e-12 * k2

  return g2_val


# Cooling due to H0 collisional excitation via collision with electrons ::: in erg.s^-1.cm^3 ::: will be multiplied by nH0 and ne later in the code.
def g3(T):

  g3_val = 7.50e-19 * np.exp(-118348/T) / (1.0 + (T/1e5)**0.5)
  
  return g3_val


# Total cooling rate in erg.s^-1.cm^-3
def Lambda(T, nH, nH0):

  nHp = nH - nH0
  ne = nHp
  
  Lamb = g1(T) * ne * nH0 # collisional ionization of hydrogen ::::::: (H0 + e ---> Hp + 2e).
       + g2(T) * ne * nHp # photo-recombination of Hp with electron :: (Hp + e ---> H0 + γ).
       + g3(T) * ne * nH0 # collisional excitaion of H0.
  
  return Lam




def n_tot(nH, nH0):
  
  nHp = nH - nH0
  ne = nHp
  
  ntot = nH0 + nHp + ne
  
  return ntot



Tgrid = np.logspace(3, 9)

res = []

for T in Tgrid:
  res.append([LambColIonH0(T, ne, nH0), k2(T)])

res = np.array(res)
a1 = res[:, 0]
a2 = res[:, 1]

print(a1/a2)

plt.scatter(np.log10(Tgrid), np.log10(a1), s = 10, color = 'k')
plt.scatter(np.log10(Tgrid), np.log10(a2), s = 10, color = 'b')
plt.show()








