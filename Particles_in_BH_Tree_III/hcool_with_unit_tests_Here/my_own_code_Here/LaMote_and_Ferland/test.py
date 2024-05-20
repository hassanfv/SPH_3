
import numpy as np
import matplotlib.pyplot as plt


# Reaction: (Hp + e ---> H0 + γ)  ::: photo-recombination of Hp and e.
def k2(T):

  k2_val = 8.4e-11 / T**0.5 / (T/1e3)**0.2 / (1.0 + (T/1e6)**0.7)

  return k2_val


# Cooling rate due to Hp photo-recombination: (Hp + e ---> H0 + γ)  ::: in erg.s^-1.cm^3  ::: will be multiplied by nHp and ne later in the code.
def g2(T):

  g2_val = 8.70e-27 * T**0.5 / (T/1e3)**0.2 / (1.0 + (T/1e6)**0.7)

  return g2_val




# Define the coefficients from the previous OCR result
a = 1.00028519708435
b = -7.569939287228937E-06
c = 2.79188868562404E-08
d = -1.289820289839189E-10
e = 7.829204293134294E-12
f = 0.2731170438382388
g = 6.086879204730784E-14
h = -0.0003748988159766978
i = 270.245476366191
j = -1.98263435534978E+09
k = -17.028197093979
l = 4.516090033327356E-05
m = 1.08832467825823

# gammaFunc gives the cooling_rate / recombination_rate ratio. Note: t = T * n**2 / Z**2 (n: the principle quantum number, Z: nuclear charge)
def gammaFunc(t):
  if t < 10**2:
    return 1.0
  
  if 10**2 <= t < 7.4 * 10**5:
    return (a + b * t + c * t**1.5 + d * t**2 + e * t**2 * np.log(t))
  
  if 7.4e5 <= t < 5e10:
    return f + g * t + h * (np.log(t))**2 + i / t**0.5 + j * np.log(t)/t**2
  
  if 5e10 <= t < 3e14:
    return 1.0 / (k + l * t**0.5 + m * np.log(t))
  
  if t >= 3e14:
    return 1.289e11 * t**-0.9705


n = 1 # Hydrogen atom
Z = 1 # Hydrogen atom

Tgrid = np.logspace(4, 9)

kB = 1.3807e-16

res = []
for T in Tgrid:

  t = T * n**2 / Z**2
  
  coeff = gammaFunc(t)
  
  CRate = coeff * k2(T) * 13.6 * 1.60219e-12
  
  res.append([T, CRate, g2(T), kB * T * k2(T)])

res = np.array(res)

T = res[:, 0]
C1 = res[:, 1]
C2 = res[:, 2]
C3 = res[:, 3]

print(np.sort(C1))

plt.scatter(np.log10(T), np.log10(C1), s = 5, color = 'k')
plt.scatter(np.log10(T), np.log10(C2), s = 5, color = 'b')
plt.scatter(np.log10(T), np.log10(C3), s = 5, color = 'r')

plt.show()



