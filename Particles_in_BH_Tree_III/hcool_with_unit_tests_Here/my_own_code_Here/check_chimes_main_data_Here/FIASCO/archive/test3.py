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

ne_min = 1e8 * u.cm**(-3)
ne_max = 1e10 * u.cm**(-3)
ne = ne_min + (ne_max - ne_min) / t_mid * t
ne[t>t_mid] = ne_max - (ne_max - ne_min) / t_mid * (t[t>t_mid] - t_mid)

temperature_array = np.logspace(4, 8, 1000) * u.K
carbon = Element('carbon', temperature_array)
func_interp = interp1d(carbon.temperature.to_value('K'), carbon.equilibrium_ionization.value,
                       axis=0, kind='cubic', fill_value='extrapolate')
carbon_ieq = u.Quantity(func_interp(Te.to_value('K')))

carbon_nei = np.zeros(t.shape + (carbon.atomic_number + 1,))
carbon_nei[0, :] = carbon_ieq[0,:]

print(carbon._rate_matrix.value.shape)
print(carbon.temperature.to_value('K').shape)


func_interp = interp1d(carbon.temperature.to_value('K'), carbon._rate_matrix.value,
                       axis=0, kind='cubic', fill_value='extrapolate')
fe_rate_matrix = func_interp(Te.to_value('K')) * carbon._rate_matrix.unit
print(fe_rate_matrix.shape)

identity = u.Quantity(np.eye(carbon.atomic_number + 1))
print(identity.shape)


for i in range(1, t.shape[0]):
    dt = t[i] - t[i-1]
    term1 = identity - ne[i] * dt/2. * fe_rate_matrix[i, ...]
    term2 = identity + ne[i-1] * dt/2. * fe_rate_matrix[i-1, ...]
    carbon_nei[i, :] = np.linalg.inv(term1) @ term2 @ carbon_nei[i-1, :]
    carbon_nei[i, :] = np.fabs(carbon_nei[i, :])
    carbon_nei[i, :] /= carbon_nei[i, :].sum()

carbon_nei = u.Quantity(carbon_nei)

print(carbon_nei.shape)

plt.figure(figsize=(12,4))
for ion in carbon:
    plt.plot(t, carbon_nei[:, ion.charge_state], ls='-', label=ion.ion_name_roman,)
plt.xlim(t[[0,-1]].value)
plt.legend(ncol=4, frameon=False)
plt.show()





