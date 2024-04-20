import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp
import pickle

# Total cooling rate in erg.s^-1.cm^-3 ===>NOTE that Hep and Hepp is excluded in free-free here as we are only considering H in this code!
def Lambda(T, nH, nH0, nHe0, nHep):

  Tx = np.log10(T)

  nHp = nH - nH0
  nHepp = nHe - nHe0 - nHep
  
  ne = nHp + nHep + 2.0 * nHepp
  
  Lamb = (
           10**g1(Tx) * ne * nH0  # H0
         + 10**g2(Tx)  * ne * nHp # Hp
         + 10**g3(Tx) * nHe0 * ne # He0 
         + 10**g4(Tx) * nHep * ne # Hep 
         + 10**g5(Tx) * nHepp * ne# Hepp
        )
  
  return Lamb


#----- n_tot
def n_tot(nH, nH0, nHe0, nHep):
  
  nHp = nH - nH0

  nHepp = nHe - nHe0 - nHep
  
  ne = nHp + nHep + 2.0 * nHepp
  
  ntot = nH0 + nHp + nHe0 + nHep + nHepp + ne
  
  return ntot



#-----------------------------------------
#------ Solution of the ODE Section ------
#-----------------------------------------
def func(t, y):

  nH0, nHe0, nHep, T = y
  
  Tx= np.log10(T)
  
  nHp = nH - nH0

  nHepp = nHe - nHe0 - nHep
  
  ne = nHp + nHep + 2.0 * nHepp
  
  dnH0_dt = 10**k2(Tx) * nHp * ne - 10**k1(Tx) * nH0 * ne
  
  dnHe0_dt = 10**k5(Tx) * nHep * ne - 10**k3(Tx) * nHe0 * ne
  
  dnHep_dt = 10**k6(Tx) * nHepp * ne + 10**k3(Tx) * nHe0 * ne - 10**k4(Tx) * nHep * ne - 10**k5(Tx) * nHep * ne
  
  Lamb = Lambda(T, nH, nH0, nHe0, nHep)
  
  ntot = n_tot(nH, nH0, nHe0, nHep)

  dT_dt = -1.0 * (gamma - 1.0) / kB / ntot * Lamb
  
  return [dnH0_dt, dnHe0_dt, dnHep_dt, dT_dt]



with h5py.File('chimes_main_data.hdf5', 'r') as file:
  rates = file['T_dependent/rates'][:]
  ratesAB = file['recombination_AB/rates_caseA'][:]
  Temp = file['TableBins/Temperatures'][:]
  cooling_rates = file['cooling/rates'][:]

# NOTE: I used "T_dependent_reactants.py" code to find indices 111, 0, 108, etc!!!

#---- Reaction rates -------
k1x = rates[111, :] # Reaction: (H0 + e ---> Hp + 2e) ::: H0 Collisional ionization
k2x = ratesAB[0, :] # Reaction: (Hp + e ---> H0 + γ)  ::: photo-recombination of Hp and e.
k3x = rates[108, :] # Reaction: (He0 + e ---> Hep + 2e) ::: He0 Collisional ionization 
k4x = rates[0, :]   # Reaction: (Hep + e ---> Hepp + 2e) ::: Hep Collisional ionization
k5x = ratesAB[1, :] # Reaction: (Hep + e ---> He0 + γ)  ::: photo-recombination of Hep and e. !!! PROBABLY di-electric is also included ?????
k6x = rates[222, :] # Reaction: (Hepp + e ---> Hep + γ)  ::: photo-recombination of Hepp and e.


#---- Cooling rates ------
g1x = cooling_rates[0, :] # cooling via H0
g2x = cooling_rates[1, :] # cooling via Hp
g3x = cooling_rates[2, :] # cooling via He0
g4x = cooling_rates[3, :] # cooling via Hep
g5x = cooling_rates[4, :] # cooling via Hepp


k1 = interp1d(Temp, k1x, kind='linear', fill_value="extrapolate")
k2 = interp1d(Temp, k2x, kind='linear', fill_value="extrapolate")
k3 = interp1d(Temp, k3x, kind='linear', fill_value="extrapolate")
k4 = interp1d(Temp, k4x, kind='linear', fill_value="extrapolate")
k5 = interp1d(Temp, k5x, kind='linear', fill_value="extrapolate")
k6 = interp1d(Temp, k6x, kind='linear', fill_value="extrapolate")

g1 = interp1d(Temp, g1x, kind='linear', fill_value="extrapolate")
g2 = interp1d(Temp, g2x, kind='linear', fill_value="extrapolate")
g3 = interp1d(Temp, g3x, kind='linear', fill_value="extrapolate")
g4 = interp1d(Temp, g4x, kind='linear', fill_value="extrapolate")
g5 = interp1d(Temp, g5x, kind='linear', fill_value="extrapolate")


gamma = 5./3.
kB = 1.3807e-16
X = 0.76
Y = 1.0 - X

nH = 1000.0

He_solar = 10**(-1.07)
nHe = He_solar * nH
print('nHe (cm^-3) = ', nHe)

y0 = [1e-3, 1.6e-6, 6e-3, 1e6]

t_span = (1*3.16e7, 20000*3.16e7)

solution = solve_ivp(func, t_span, y0, method='LSODA', dense_output=True)

t = np.linspace(t_span[0], t_span[1], 10000) # This 10000 is not years, it is the number of points in linspace !!!!
y = solution.sol(t)

print(y.shape)

t_yrs = t / 3.16e7

nH0 = y[0, :]
nHe0 = y[1, :]
nHep = y[2, :]

nHepp = nHe - nHe0 - nHep

T = y[3, :]

nHp = (nH - nH0)

print('tot He = ', nHe0 + nHep + nHepp)

#----- Preparing cooling rate for plotting -----
res = []
for Tx, nH0x, nHe0x, nHepx in zip(T, nH0, nHe0, nHep):

  lmb = Lambda(Tx, nH, nH0x, nHe0x, nHepx)
  
  res.append([Tx, lmb])

res = np.array(res)

Tx = res[:, 0]
lmb = res[:, 1]


#------ Result from "test_primordial_hdf5_v2.py" code -----
with open('chimesRes.pkl', 'rb') as f:
  df = pickle.load(f)
# dictx = {'t_Arr_in_yrs': t_Arr_in_yrs, 'TEvol': TEvol, 'nHe0': nHe0, 'nHep': nHep, 'nHepp': nHepp}
t_Arr_in_yrsx = df['t_Arr_in_yrs']
TEvolx = df['TEvol']
nH0x = df['nH0']
nHpx = df['nHp']
nHe0x = df['nHe0']
nHepx = df['nHep']
nHeppx = df['nHepp']
nHeTotx = nHe0x + nHepx + nHeppx
#----------------------------------------------------------


plt.figure(figsize = (16, 8))

plt.subplot(2, 3, 1)
plt.scatter(t_yrs, np.log10(T), s = 5, color = 'k', label = 'my own code')
plt.scatter(t_Arr_in_yrsx, np.log10(TEvolx), s = 2, color = 'orange', label = 'chimes result', linestyle = '--')
plt.xlim(0, 3000)
plt.ylim(3, 8)
plt.legend()

plt.subplot(2, 3, 2)

plt.plot(t_yrs, nHe0/nHe, color = 'r', label = 'nHe0')
plt.plot(t_yrs, nHep/nHe, color = 'g', label = 'nHep')
plt.plot(t_yrs, nHepp/nHe, color = 'b', label = 'nHepp')

plt.plot(t_Arr_in_yrsx, nHe0x/nHeTotx, color = 'r', label = 'nHe0 - chimes', linestyle = ':')
plt.plot(t_Arr_in_yrsx, nHepx/nHeTotx, color = 'g', label = 'nHep - chimes', linestyle = ':')
plt.plot(t_Arr_in_yrsx, nHeppx/nHeTotx, color = 'b', label = 'nHepp - chimes', linestyle = ':')

plt.xlim(0, 5000)
plt.ylim(1e-7, 2.)

plt.yscale('log')
plt.title('solve_ivp')
plt.legend()


plt.subplot(2, 3, 3)
nHeTot = nHe0 + nHep + nHepp
plt.plot(T, nH0/nH, label = 'nH0', color = 'r')
plt.plot(T, nHp/nH, label = 'nHp', color = 'g')
plt.plot(T, nHe0/nHeTot, label = 'nHe0', color = 'b')
plt.plot(T, nHep/nHeTot, label = 'nHep', color = 'orange')
plt.plot(T, nHepp/nHeTot,label = 'nHepp', color = 'purple')

plt.plot(TEvolx, nH0x/nH, label = 'nH0 - chimes', color = 'r', linestyle = ':')
plt.plot(TEvolx, nHpx/nH, label = 'nHp - chimes', color = 'g', linestyle = ':')
plt.plot(TEvolx, nHe0x/nHeTotx, label = 'nHe0 - chimes', color = 'b', linestyle = ':')
plt.plot(TEvolx, nHepx/nHeTotx, label = 'nHep - chimes', color = 'orange', linestyle = ':')
plt.plot(TEvolx, nHeppx/nHeTotx,label = 'nHepp - chimes', color = 'purple', linestyle = ':')

plt.yscale('log')
plt.xscale('log')
plt.ylim(2e-3, 1.2)
plt.xlim(1e4, 5e6)
plt.legend()

plt.subplot(2, 3, 4)
plt.scatter(np.log10(Tx), np.log10(lmb/nH/nH), s = 5, color = 'k')
plt.xlim(3.5, 8.25)
plt.ylim(-25, -21.5)

plt.tight_layout()

plt.savefig('HeH_hfv.png')

plt.show()




