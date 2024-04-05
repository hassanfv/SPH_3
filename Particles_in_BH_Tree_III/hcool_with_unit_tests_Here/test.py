
from photolibs3 import *
import numpy as np
import matplotlib.pyplot as plt


XH = 0.7
kB = 1.3807e-16  # cm2 g s-2 K-1
mH = 1.6726e-24 # gram
gamma = 5.0/3.0

mu = 0.6

T = 1e5 # K
u_old = kB * T / (gamma - 1.0) / mH / mu

print(f'u_old = {u_old:.5E}')

nH = 100.0

rho = nH * mH / XH

print(f'rho = {rho:.4E}')

dt_yrs = 10.0

dt_sec = dt_yrs * 365.25 * 24.0 * 3600

print(f'dt_sec = {dt_sec:.3E}')

u_after = DoCooling_h(rho, u_old, dt_sec, XH)


print(f'u_after = {u_after:.5E}')
print()
print(f'T before = {T}')
print(f'T after = {(u_after/kB * (gamma - 1.0) * mH * mu):.2f}')


res = []
for dt in np.arange(1, dt_yrs+1):

  dt_s = dt * 365.25 * 24.0 * 3600

  u_after = DoCooling_h(rho, u_old, dt_s, XH)

  T1 = u_after/kB * (gamma - 1.0) * mH * mu
  res.append(T1)
  
  print(dt)

res = np.array(res)

plt.plot(np.arange(len(res)), res)
plt.xlim(0, 10)
plt.ylim(9e4, 1e5)
plt.show()






