
import numpy as np
import matplotlib.pyplot as plt


# Reaction: (Hp + e ---> H0 + γ)  ::: photo-recombination of Hp and e.
def k2(T):

  k2_val = 8.4e-11 / T**0.5 / (T/1e3)**0.2 / (1.0 + (T/1e6)**0.7)

  return k2_val
#--------------


# Cooling rate due to Hp photo-recombination: (Hp + e ---> H0 + γ)  ::: in erg.s^-1.cm^3  ::: will be multiplied by nHp and ne later in the code.
def g2(T):

  g2_val = 8.70e-27 * T**0.5 / (T/1e3)**0.2 / (1.0 + (T/1e6)**0.7)

  return g2_val


#----- gx
def gx(T):

  gx_val = 13.5984 * 1.60219e-12 * k2(T)

  return gx_val



Tgrid = np.logspace(4, 8)

res = []
for T in Tgrid:

  res.append([g2(T), gx(T)])
  
res = np.array(res)
L1 = res[:, 0]
Lx = res[:, 1]

plt.scatter(np.log10(Tgrid), np.log10(L1), s = 5, color = 'k')
plt.scatter(np.log10(Tgrid), np.log10(Lx), s = 5, color = 'b')

plt.show()




