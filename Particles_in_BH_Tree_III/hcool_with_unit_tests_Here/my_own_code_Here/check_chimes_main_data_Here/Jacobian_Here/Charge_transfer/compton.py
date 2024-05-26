
import numpy as np
import matplotlib.pyplot as plt


ne = 100.0
z = 0.0


T = np.logspace(4, 8)

L1 = 5.41e-36 * ne * T * (1.0 + z)**4 # Compton cooling - Katz famous paper!

TCMB_0 = 2.7255
TCMB = TCMB_0 * (1.0 + z)
L2 = 1.017e-37 * TCMB**4 * (T - TCMB) * ne

plt.loglog(T, L1, label = 'Katz')
plt.loglog(T, L2, label = 'new')

plt.legend()

plt.show()

