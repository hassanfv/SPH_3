import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp
from scipy.interpolate import RegularGridInterpolator
import pickle

# Total cooling rate in erg.s^-1.cm^-3 ===>NOTE that Hep and Hepp is excluded in free-free here as we are only considering H in this code!
def Lambda(T, nH, nH0, nHe0, nHep, nC0, nC1, nC2, nC3, nC4, nC5, nC6): #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Tx = np.log10(T)

  nHp = nH - nH0
  nHepp = nHe - nHe0 - nHep
  
  ne = nHp + nHep + 2.0 * nHepp + (nC1 + 2.0 * nC2 + 3. * nC3 + 4. * nC4 + 5. * nC5 + 6. * nC6) #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  Lamb = (
           10**g1(Tx) * ne * nH0  # H0
         + 10**g2(Tx)  * ne * nHp # Hp
         + 10**g3(Tx) * nHe0 * ne # He0 
         + 10**g4(Tx) * nHep * ne # Hep 
         + 10**g5(Tx) * nHepp * ne# Hepp
         + 10**C0_cooling_rate(T, nH0, ne, nHp, Temp_4d, HIDensity_4d, elecDensity_4d, HIIDensity_4d) * nC0 * ne # cooling via C0
         + 10**Cp_cooling_rate(T, ne, Temp_2d, elecDensity_2d) * nC1 * ne # cooling via Cp or C1
         + 10**gC2(Tx) * nC2 * ne # C2
         + 10**gC3(Tx) * nC3 * ne # C3
         + 10**gC4(Tx) * nC4 * ne # C4
         + 10**gC5(Tx) * nC5 * ne # C5
         + 10**gC6(Tx) * nC6 * ne # C6
        )
  
  return Lamb


#----- n_tot
def n_tot(nH, nH0, nHe0, nHep, nC, nC0, nC1, nC2, nC3, nC4, nC5, nC6): #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  nHp = nH - nH0

  nHepp = nHe - nHe0 - nHep
  
  ne = nHp + nHep + 2.0 * nHepp + (nC1 + 2. * nC2 + 3. * nC3 + 4. * nC4 + 5. * nC5 + 6. * nC6) #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ntot = nH0 + nHp + nHe0 + nHep + nHepp + ne + nC
  
  return ntot



#-----------------------------------------
#------ Solution of the ODE Section ------
#-----------------------------------------
def func(t, y):

  nH0, nHe0, nHep, nC0, nC1, nC2, nC3, nC4, nC5, nC6, T = y #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  Tx= np.log10(T)
  
  nHp = nH - nH0

  nHepp = nHe - nHe0 - nHep
  
  #nC3 = nC - (nC0 + nC1 + nC2)
  
  ne = nHp + nHep + 2.0 * nHepp + (nC1 + 2. * nC2 + 3. * nC3 + 4. * nC4 + 5. * nC5 + 6. * nC6) #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  dnH0_dt = 10**k2(Tx) * nHp * ne - 10**k1(Tx) * nH0 * ne
  
  dnHe0_dt = 10**k5(Tx) * nHep * ne - 10**k3(Tx) * nHe0 * ne
  
  dnHep_dt = 10**k6(Tx) * nHepp * ne + 10**k3(Tx) * nHe0 * ne - 10**k4(Tx) * nHep * ne - 10**k5(Tx) * nHep * ne
  
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  dnC0_dt = 10**C_1_0(Tx) * ne * nC1 - 10**C_0_1(Tx) * ne * nC0
  dnC1_dt = 10**C_0_1(Tx) * ne * nC0 + 10**C_2_1(Tx) * ne * nC2 - 10**C_1_0(Tx) * ne * nC1 - 10**C_1_2(Tx) * ne * nC1
  dnC2_dt = 10**C_1_2(Tx) * ne * nC1 + 10**C_3_2(Tx) * ne * nC3 - 10**C_2_1(Tx) * ne * nC2 - 10**C_2_3(Tx) * ne * nC2
  dnC3_dt = 10**C_2_3(Tx) * ne * nC2 + 10**C_4_3(Tx) * ne * nC4 - 10**C_3_2(Tx) * ne * nC3 - 10**C_3_4(Tx) * ne * nC3
  dnC4_dt = 10**C_3_4(Tx) * ne * nC3 + 10**C_5_4(Tx) * ne * nC5 - 10**C_4_3(Tx) * ne * nC4 - 10**C_4_5(Tx) * ne * nC4
  dnC5_dt = 10**C_4_5(Tx) * ne * nC4 + 10**C_6_5(Tx) * ne * nC6 - 10**C_5_4(Tx) * ne * nC5 - 10**C_5_6(Tx) * ne * nC5
  dnC6_dt = 10**C_5_6(Tx) * ne * nC5 - 10**C_6_5(Tx) * ne * nC6
  
  Lamb = Lambda(T, nH, nH0, nHe0, nHep, nC0, nC1, nC2, nC3, nC4, nC5, nC6) #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ntot = n_tot(nH, nH0, nHe0, nHep, nC, nC0, nC1, nC2, nC3, nC4, nC5, nC6) #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  dT_dt = -1.0 * (gamma - 1.0) / kB / ntot * Lamb
  
  return [dnH0_dt, dnHe0_dt, dnHep_dt, dnC0_dt, dnC1_dt, dnC2_dt, dnC3_dt, dnC4_dt, dnC5_dt, dnC6_dt, dT_dt] #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



#----- C0_cooling_rate 
def C0_cooling_rate(T, nHI, nelec, nHII, Temp_4d, HIDensity_4d, elecDensity_4d, HIIDensity_4d):

  T = np.log10(T)
  nHI = np.log10(nHI)
  nelec = np.log10(nelec)
  nHII = np.log10(nHII)

  if T <= 4:
    C0_rates = rates_4d[0, :]
    interp_4d = RegularGridInterpolator((Temp_4d, HIDensity_4d, elecDensity_4d, HIIDensity_4d), C0_rates)
    res = -25.0#interp_4d(np.array([T, nHI, nelec, nHII]))[0]
  else:
    C0_rates = rates_hiT_4d[0, :]
    interp_4d = interp1d(Temp_hiT_4d, C0_rates, kind='linear', fill_value="extrapolate")
    res = interp_4d(T)

  return res


#----- Cp_cooling_rate 
def Cp_cooling_rate(T, nelec, Temp_2d, elecDensity_2d): # include Temp_hiT here !!!!!!!!!!!

  T = np.log10(T)
  nelec = np.log10(nelec)

  if T <= 4:
    Cp_rates = rates_2d[0, :]
    interp_2d = RegularGridInterpolator((Temp_2d, elecDensity_2d), Cp_rates)
    res = interp_2d(np.array([T, nelec]))[0]
  else:
    Cp_rates = rates_hiT_2d[0, :]
    interp_2d = interp1d(Temp_hiT_2d, Cp_rates, kind='linear', fill_value="extrapolate")
    res = interp_2d(T)

  return res


with h5py.File('chimes_main_data.hdf5', 'r') as file:
  rates = file['T_dependent/rates'][:]
  ratesAB = file['recombination_AB/rates_caseA'][:]
  Temp = file['TableBins/Temperatures'][:]
  cooling_rates = file['cooling/rates'][:]
  
  #------------------------------------------------------------
  #------------------- Carbon Section -------------------------
  #------------------------------------------------------------
  Temp_4d = file['TableBins/cool_4d_Temperatures'][:]
  HIDensity_4d = file['TableBins/cool_4d_HIDensities'][:] # only used for low T
  elecDensity_4d = file['TableBins/cool_4d_ElectronDensities'][:] # only used for low T
  HIIDensity_4d = file['TableBins/cool_4d_HIIDensities'][:] # only used for low T
  rates_4d = file['cooling/rates_4d'][:] # NOTE it is rates_4d for low Temp
  #-------- hiT_4d ---------------
  Temp_hiT_4d = file['TableBins/cool_hiT_4d_Temperatures'][:]
  rates_hiT_4d = file['cooling/rates_hiT_4d'][:] # NOTE it is rates_4d for high Temp
  #------------------ CII section --------------------------
  coolants_2d = file['cooling/coolants_2d'][:] # is the same for low and hiT.... Only Temp will change!! And in hiT it is only a function of T!!!
  Temp_2d = file['TableBins/cool_2d_Temperatures'][:]
  elecDensity_2d = file['TableBins/cool_2d_ElectronDensities'][:]
  rates_2d = file['cooling/rates_2d'][:] # NOTE it is rates_2d for low Temp
  #-------- hiT_2d ---------------
  Temp_hiT_2d = file['TableBins/cool_hiT_2d_Temperatures'][:]
  rates_hiT_2d = file['cooling/rates_hiT_2d'][:] # NOTE it is rates_2d for high Temp



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



C_0_1x = rates[220, :]
C_1_0x = rates[216, :]
C_1_2x = rates[217, :]
C_2_1x = rates[212, :]
C_2_3x = rates[213, :]
C_3_2x = rates[209, :]
C_3_4x = rates[210, :]
C_4_3x = rates[223, :]
C_4_5x = rates[206, :]
C_5_4x = rates[243, :]
C_5_6x = rates[226, :]
C_6_5x = rates[241, :]


plt.plot(Temp, C_1_2x, color = 'lime')
plt.ylim(-12, -6)
plt.show()
s()


C_0_1 = interp1d(Temp, C_0_1x, kind='linear', fill_value="extrapolate")
C_1_0 = interp1d(Temp, C_1_0x, kind='linear', fill_value="extrapolate")
C_1_2 = interp1d(Temp, C_1_2x, kind='linear', fill_value="extrapolate")
C_2_1 = interp1d(Temp, C_2_1x, kind='linear', fill_value="extrapolate")
C_2_3 = interp1d(Temp, C_2_3x, kind='linear', fill_value="extrapolate")
C_3_2 = interp1d(Temp, C_3_2x, kind='linear', fill_value="extrapolate")
C_3_4 = interp1d(Temp, C_3_4x, kind='linear', fill_value="extrapolate")
C_4_3 = interp1d(Temp, C_4_3x, kind='linear', fill_value="extrapolate")
C_4_5 = interp1d(Temp, C_4_5x, kind='linear', fill_value="extrapolate")
C_5_4 = interp1d(Temp, C_5_4x, kind='linear', fill_value="extrapolate")
C_5_6 = interp1d(Temp, C_5_6x, kind='linear', fill_value="extrapolate")
C_6_5 = interp1d(Temp, C_6_5x, kind='linear', fill_value="extrapolate")

# cooling via C0 and C1 are determined using "C0_cooling_rate" and "C0_cooling_rate" functions inside the "Lambda" function directly!!
gC2x = cooling_rates[5, :]
gC3x = cooling_rates[6, :]
gC4x = cooling_rates[7, :]
gC5x = cooling_rates[8, :]
gC6x = cooling_rates[9, :]


gC2 = interp1d(Temp, gC2x, kind='linear', fill_value="extrapolate")
gC3 = interp1d(Temp, gC3x, kind='linear', fill_value="extrapolate")
gC4 = interp1d(Temp, gC4x, kind='linear', fill_value="extrapolate")
gC5 = interp1d(Temp, gC5x, kind='linear', fill_value="extrapolate")
gC6 = interp1d(Temp, gC6x, kind='linear', fill_value="extrapolate")




gamma = 5./3.
kB = 1.3807e-16
X = 0.76
Y = 1.0 - X

nH = 1000.0

C_solar = 10**(-3.57)
nC = C_solar * nH

He_solar = 10**(-1.07)
nHe = He_solar * nH
print('nHe (cm^-3) = ', nHe)

#y0 = [1e-6, 0.6e-6, 0e-3, 0e-7, 0e-7, 0e-4, 0e-4, 0e-2, 0e-2, 0e-2, 1e6] #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
y0 = [1e-2, 1.6e-3, 6e-2, nC/8, 0/8, 0/8, 0/8, 0/8, 0/8, 0/8, 1e6] #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

t_span = (1*3.16e7, 3000*3.16e7)

solution = solve_ivp(func, t_span, y0, method='LSODA', dense_output=True)

t = np.linspace(t_span[0], t_span[1], 3000) # This 10000 is not years, it is the number of points in linspace !!!!
y = solution.sol(t)

print(y.shape)

t_yrs = t / 3.16e7

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
nH0 = y[0, :]
nHe0 = y[1, :]
nHep = y[2, :]
nC0 = y[3, :]
nC1 = y[4, :]
nC2 = y[5, :]
nC3 = y[6, :]
nC4 = y[7, :]
nC5 = y[8, :]
nC6 = y[9, :]
T = y[10, :]

print(nC0)
print()
print(nC1)
print()
print(nC2)
print()
print(nC3)
print()
print(nC4)
print()
print(nC5)
print()
print(nC6)
print()
print('orig nC = ', nC)
print()
print('tot_nC_after = ', nC0+nC1+nC2+nC3+nC4+nC5+nC6)
print()
print(T)
print()
#s()

nHepp = nHe - nHe0 - nHep


nHp = (nH - nH0)

print('orig nH = ', nH)

#----- Preparing cooling rate for plotting -----
res = []
for Tx, nH0x, nHe0x, nHepx, nC0x, nC1x, nC2x, nC3x, nC4x, nC5x, nC6x in zip(T, nH0, nHe0, nHep, nC0, nC1, nC2, nC3, nC4, nC5, nC6): #!!!!!!!!!!!!!!!!!!!!!

  lmb = Lambda(Tx, nH, nH0x, nHe0x, nHepx, nC0x, nC1x, nC2x, nC3x, nC4x, nC5x, nC6x) #!!!!!!!!!!!!!!!!!!!!!!!!!!
  #print(Tx, nH, nH0x, nHe0x, nHepx, nC0x, nC1x, nC2x, nC3x, nC4x, nC5x, nC6x, lmb)
  
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
nHx = nH0x + nHpx
nHe0x = df['nHe0']
nHepx = df['nHep']
nHeppx = df['nHepp']
nHeTotx = nHe0x + nHepx + nHeppx

nC0x = df['nC0']
nC1x = df['nC1']
nC2x = df['nC2']
nC3x = df['nC3']
nC4x = df['nC4']
nC5x = df['nC5']
nC6x = df['nC6']
nCx = nC0x + nC1x + nC2x + nC3x + nC4x + nC5x + nC6x
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
plt.xlim(1e4, 1e6)
plt.legend()

plt.subplot(2, 3, 4)
plt.scatter(np.log10(Tx), np.log10(lmb/nH/nH), s = 5, color = 'k')
#plt.xlim(3.5, 8.25)
plt.ylim(-25, -21.5)

plt.subplot(2, 3, 5)
plt.plot(T, nC0/nC, label = 'nC0', color = 'r')
plt.plot(T, nC1/nC, label = 'nC1', color = 'g')
plt.plot(T, nC2/nC, label = 'nC2', color = 'b')
plt.plot(T, nC3/nC, label = 'nC3', color = 'orange')
plt.plot(T, nC4/nC, label = 'nC4', color = 'purple')
plt.plot(T, nC5/nC, label = 'nC5', color = 'lime')
plt.plot(T, nC6/nC, label = 'nC6', color = 'pink')

plt.plot(TEvolx, nC0x/nCx, label = 'nC0', color = 'r', linestyle = ':')
plt.plot(TEvolx, nC1x/nCx, label = 'nC1', color = 'g', linestyle = ':')
plt.plot(TEvolx, nC2x/nCx, label = 'nC2', color = 'b', linestyle = ':')
plt.plot(TEvolx, nC3x/nCx, label = 'nC3', color = 'orange', linestyle = ':')
plt.plot(TEvolx, nC4x/nCx, label = 'nC4', color = 'purple', linestyle = ':')
plt.plot(TEvolx, nC5x/nCx, label = 'nC5', color = 'lime', linestyle = ':')
plt.plot(TEvolx, nC6x/nCx, label = 'nC6', color = 'pink', linestyle = ':')

plt.yscale('log')
plt.xscale('log')
plt.ylim(2e-3, 1.2)
plt.xlim(1e4, 1e6)
plt.legend()

plt.tight_layout()

plt.savefig('HeH_hfv.png')

plt.show()




