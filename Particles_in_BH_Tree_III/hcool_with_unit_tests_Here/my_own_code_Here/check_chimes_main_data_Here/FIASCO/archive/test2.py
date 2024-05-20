import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.visualization import quantity_support
from scipy.interpolate import interp1d
from fiasco import Element



t = np.arange(0, 30, 0.01) * u.s
Te_min = 1e4 * u.K
Te_max = 2e6 * u.K
t_mid = 15 * u.s
Te = Te_min + (Te_max - Te_min) / t_mid * t
Te[t>t_mid] = Te_max - (Te_max - Te_min) / t_mid * (t[t>t_mid] - t_mid)

temperature_array = np.logspace(4, 8, 1000) * u.K
carbon = Element('carbon', temperature_array)
func_interp = interp1d(carbon.temperature.to_value('K'), carbon.equilibrium_ionization.value,
                       axis=0, kind='cubic', fill_value='extrapolate')
carbon_ieq = u.Quantity(func_interp(Te.to_value('K')))



for ion in carbon:
    plt.plot(Te, carbon_ieq[:, ion.charge_state], label=ion.ion_name_roman)
plt.legend(ncol=4, frameon=False)
plt.xscale('log')
plt.yscale('log')
plt.ylim(10**-4.5, 2)
plt.show()
