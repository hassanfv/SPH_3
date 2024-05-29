
import numpy as np


gamma = 5./3.
kB = 1.3807e-16


T = 1e6 # K

nH = 1000.0

He_solar = 10**(-1.0) # See Table_1 in Wiersma et al - 2009, 393, 99–107
nHe = He_solar * nH
print('nHe (cm^-3) = ', nHe)

C_solar = 10**(-3.61) # See Table_1 in Wiersma et al - 2009, 393, 99–107
nC = C_solar * nH

print('nC (cm^-3) = ', nC)

ntot = 1200 # educated guess

u = 1. / (gamma - 1.) / ntot / kB * T

print(f'u = {u:.3E}')


