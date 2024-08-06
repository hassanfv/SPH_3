
import numpy as np


r = np.arange(0.01, 0.71, 0.1)
L = np.arange(12., 23., 1.)
t = np.arange(0, 10000., 1.)

Nr = len(r)
NL = len(L)
Nt = len(t)

T = np.zeros((Nr, NL, Nt))

for i in range(Nr):
  for j in range(NL):
    for k in range(Nt):
      
      T[i, j, k] = 10. + (Nt - k) * j * i

print(T)

