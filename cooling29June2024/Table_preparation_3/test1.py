
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


df = pd.read_csv('nH_vs_time.csv')
# ['nH', 'n_iter', 'dt']
print(df.columns)

nH = df['nH'].values
n_iter = df['n_iter'].values
dt = df['dt'].values

t = n_iter * dt
t_yrs = t / 3600. / 24. / 365.25

#--- Interpolation ------
nHG = np.linspace(-4.0, 4.0, 81)  # Generates 81 points from -4 to 4 inclusively
f_linear = interp1d(nH, np.log10(t_yrs), bounds_error=False, fill_value="extrapolate")
log_tG = f_linear(nHG)

jj = 65
nH_i, log_t_i = nHG[jj], log_tG[jj]
print(nH_i, log_t_i)

t_i = 10**log_t_i * 3600. * 24. * 365.25 # converting back to seconds!
dt_i = t_i / 30000
print(f'dt_i = {dt_i:.2E}')
#------------------------

plt.scatter(nH, np.log10(t_yrs), color = 'k')
plt.scatter(nHG, log_tG, color = 'b', s = 10)

plt.show()


