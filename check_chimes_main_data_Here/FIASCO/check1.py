import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.visualization import quantity_support
from fiasco import Ion



Te = np.geomspace(1e4, 1e8, 101) * u.K
ne = 1e8 * u.cm**-3

ion = Ion('H I', Te)
print(ion)


contribution_func = ion.contribution_function(ne)
wlen = 1031.93 * u.Angstrom

print('Te.shape = ', Te.shape)
print(contribution_func.shape)


transitions = ion.transitions.wavelength[~ion.transitions.is_twophoton]
idx = np.argmin(np.abs(transitions - wlen))


Lamb = np.sum(contribution_func, axis=2)


plt.plot(Te, contribution_func[:, 0, idx], label=f'{ion.atomic_symbol} {ion.charge_state}+ {wlen}')

plt.plot(Te, Lamb, color = 'r')

plt.title('Contribution function')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.show()



