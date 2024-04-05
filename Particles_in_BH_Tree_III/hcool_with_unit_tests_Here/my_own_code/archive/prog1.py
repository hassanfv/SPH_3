
import numpy as np
import matplotlib.pyplot as plt


#===== Gam_e_H0
def Gam_e_H0(T):
  return 5.85e-11 * T**0.5 * np.exp(-157809.1/T) / (1.0 + (T/1e5)**0.5)


# Reaction: H0 + e = Hp + 2e
def k1(T):

    T = T * 8.61732814974493E-05 # convert Kelvin to eV.
    ln_T = np.log(T)
    
    k1 = np.exp(-32.71396786
                 + 13.5365560 * ln_T
                 - 5.73932875 * ln_T**2
                 + 1.56315498 * ln_T**3
                 - 0.287705600 * ln_T**4
                 + 3.48255977 * 10**-2 * ln_T**5
                 - 2.63197617 * 10**-3 * ln_T**6
                 + 1.11954395 * 10**-4 * ln_T**7
                 - 2.03914985 * 10**-6 * ln_T**8)
    return k1


Tgrid = np.logspace(3, 9)

res = []

for T in Tgrid:
  res.append([Gam_e_H0(T), k1(T)])

res = np.array(res)
G1 = res[:, 0]
G2 = res[:, 1]

plt.scatter(np.log10(Tgrid), np.log10(G1), s = 10, color = 'k')
plt.scatter(np.log10(Tgrid), np.log10(G2), s = 10, color = 'b')
plt.show()






