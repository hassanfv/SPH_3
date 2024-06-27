import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.visualization import quantity_support
from fiasco import Ion



Te = np.geomspace(1e4, 1e8, 101) * u.K
ne = 1e1 * u.cm**-3

ion = Ion('H I', Te)
contribution_func = ion.contribution_function(ne)
Lamb_H = np.sum(contribution_func, axis=2)
Lamb_H = Lamb_H.reshape(Lamb_H.shape[0],)
Lamb_H = Lamb_H.value


ion2 = Ion('He II', Te)
contribution_func = ion2.contribution_function(ne)
Lamb_He = np.sum(contribution_func, axis=2)
Lamb_He = Lamb_He.reshape(Lamb_H.shape[0],)
Lamb_He = Lamb_He.value


Te = Te.value


plt.plot(np.log10(Te), np.log10(Lamb_H + Lamb_He), color = 'b')

plt.title('Contribution function')
#plt.xscale('log')
#plt.yscale('log')

plt.xlim(4, 8)
plt.ylim(-25, -21)

plt.legend()
plt.show()



