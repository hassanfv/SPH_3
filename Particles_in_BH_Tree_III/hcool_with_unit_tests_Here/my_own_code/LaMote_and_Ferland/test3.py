
import numpy as np
import matplotlib.pyplot as plt


# Reaction: (He0 + e ---> Hep + 2e) ::: He0 Collisional ionization 
def k3(T):

  k3_val = 2.38e-11 * T**0.5 * np.exp(-285335.4/T) / (1.0 + (T/1e5)**0.5)
  
  return k3_val
#--------------


# Cooling due to He0 collisional ioniozation ::: will be multiplied by nHe0 and ne later in the code
def g6(T):
  
  g6_val = 9.38e-22 * T**0.5 * np.exp(-285335.4/T) / (1. + (T/1e5)**0.5)
  
  return g6_val


IP_Hep = 24.59 * 1.60218e-12 # erg



Tgrid = np.logspace(4, 9)

kB = 1.3807e-16

res = []
for T in Tgrid:
  
  res.append([T, g6(T), IP_Hep * k3(T)])

res = np.array(res)

T = res[:, 0]
C1 = res[:, 1]
C2 = res[:, 2]

print(np.sort(C1))

plt.scatter(np.log10(T), np.log10(C1), s = 20, color = 'b')
plt.scatter(np.log10(T), np.log10(C2), s = 5, color = 'r')

plt.show()



