import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.visualization import quantity_support
from fiasco import Ion



ion = Ion('He III', np.logspace(4, 8, 100) * u.K)

print(ion)


attributes_and_methods = dir(ion)
print(attributes_and_methods)


plt.plot(ion.temperature, ion.recombination_rate, label='Recombination', color='C0',) # Total
plt.plot(ion.temperature, ion.ionization_rate, label='Ionization', color='C1') # Total

plt.xscale('log')
plt.yscale('log')
plt.xlim(1e4, 1e8)
plt.ylim(1e-18, 1e-6)
plt.legend()
plt.show()




