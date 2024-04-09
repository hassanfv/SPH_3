
import numpy as np
import matplotlib.pyplot as plt


# Reaction: (Hep + e ---> He0 + γ)  ::: photo-recombination of Hep and e.
def k5(T):

  k5_val = 1.5e-10 / T**0.6353
  
  return k5_val
#--------------


# Cooling via Hep recombination with electron ::: will be multiplied by nHep and ne later in the code
def g8(T):

  g8_val = 1.55e-26 * T**0.3647
  
  return g8_val



Tgrid = np.logspace(4, 9)

kB = 1.3807e-16

res = []
for T in Tgrid:
  
  res.append([T, g8(T), kB * T * k5(T)])

res = np.array(res)

T = res[:, 0]
C1 = res[:, 1]
C2 = res[:, 2]

print(np.sort(C1))

plt.scatter(np.log10(T), np.log10(C1), s = 5, color = 'k')
plt.scatter(np.log10(T), np.log10(C2), s = 5, color = 'b')

plt.show()



