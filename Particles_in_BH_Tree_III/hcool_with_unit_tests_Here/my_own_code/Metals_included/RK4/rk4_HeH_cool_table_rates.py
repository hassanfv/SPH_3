
#Ref: https://www.youtube.com/watch?v=vNoFdtcPFdk

import numpy as np
import matplotlib.pyplot as plt
import pickle

#----------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
#------------------------------------- Reaction Rates (cm^3.s^-1) -----------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
# Reaction: (H0 + e ---> Hp + 2e) ::: H0 Collisional ionization
def k1(T):
  
  k1_val = 5.85e-11 * T**0.5 * np.exp(-157809.1/T) / (1.0 + (T/1e5)**0.5)

  return k1_val
#--------------

# Reaction: (Hp + e ---> H0 + γ)  ::: photo-recombination of Hp and e.
def k2(T):

  k2_val = 8.4e-11 / T**0.5 / (T/1e3)**0.2 / (1.0 + (T/1e6)**0.7)

  return k2_val
#--------------

# Reaction: (He0 + e ---> Hep + 2e) ::: He0 Collisional ionization 
def k3(T):

  k3_val = 2.38e-11 * T**0.5 * np.exp(-285335.4/T) / (1.0 + (T/1e5)**0.5)
  
  return k3_val
#--------------

# Reaction: (Hep + e ---> Hepp + 2e) ::: Hep Collisional ionization 
def k4(T):

  k4_val = 5.68e-12 * T**0.5 * np.exp(-631515./T) / (1.0 + (T/1e5)**0.5)
  
  return k4_val
#--------------

# Reaction: (Hep + e ---> He0 + γ)  ::: photo-recombination of Hep and e.
def k5(T):

  k5_val = 1.5e-10 / T**0.6353
  
  return k5_val
#--------------

# Reaction: (Hepp + e ---> Hep + γ)  ::: photo-recombination of Hepp and e.
def k6(T):

  k6_val = 3.36e-10 / T**0.5 / (T/1e3)**0.2 / (1. + (T/1e6)**0.7)
  
  return k6_val
#--------------

#Reaction: Hep di-electric recombination (Hep + e ---> He0 + γ)
def k7(T):

  k7_val = 1.9e-3 / T**1.5 * np.exp(-470000./T) * (1.0 + 0.3 * np.exp(-94000./T))
  
  return k7_val
#--------------


# Total cooling rate in erg.s^-1.cm^-3 ===>NOTE that Hep and Hepp is excluded in free-free here as we are only considering H in this code!
def Lambda(T, nH, nH0, nHe0, nHep):

  nHp = nH - nH0
  nHepp = nHe - nHe0 - nHep
  
  ne = nHp + nHep + 2.0 * nHepp
  
  Lamb = (
           g_interp(Tref, T, g1) * ne * nH0 # collisional ionization of hydrogen ::::::: (H0 + e ---> Hp + 2e).
         + g_interp(Tref, T, g2)  * ne * nHp # photo-recombination of Hp with electron :: (Hp + e ---> H0 + γ).
         + g_interp(Tref, T, g3)  * ne * nH0 # collisional excitaion of H0.
         + g_interp(Tref, T, g4) * ne * (nHp + nHep + 4.0 * nHepp) # free-free emission
         + g_interp(Tref, T, g5) * nHep * ne # collisional excitation of Hep.
         + g_interp(Tref, T, g6) * nHe0 * ne # He0 collisional ionization
         + g_interp(Tref, T, g7) * nHep * ne # Hep collisional ionization
         + g_interp(Tref, T, g8) * nHep * ne # Hep recombination to He0
         + g_interp(Tref, T, g9) * nHepp * ne# Hepp recombination to Hep
         + g_interp(Tref, T, g10) * nHep * ne# Hep di-electric recombination to He0
        )
  
  return Lamb



#----- n_tot
def n_tot(nH, nH0, nHe0, nHep):
  
  nHp = nH - nH0

  nHepp = nHe - nHe0 - nHep
  
  ne = nHp + nHep + 2.0 * nHepp
  
  ntot = nH0 + nHp + nHe0 + nHep + nHepp + ne
  
  return ntot



#===== rk4SingleStep
def rk4SingleStep(func, dt, t0, y0):

  f1 = func(t0, y0)
  f2 = func(t0 + dt/2.0, y0 + (dt/2.0) * f1)
  f3 = func(t0 + dt/2.0, y0 + (dt/2.0) * f2)
  f4 = func(t0 + dt, y0 + dt * f3)
  
  yout = y0 + (dt / 6.0) * (f1 + 2.0 * f2 + 2.0 * f3 + f4)
  
  #print(yout)
  
  return yout



#-----------------------------------------
#------ Solution of the ODE Section ------
#-----------------------------------------
def func(t, y):

  nH0, nHe0, nHep, T = y
  
  nHp = nH - nH0

  nHepp = nHe - nHe0 - nHep
  
  ne = nHp + nHep + 2.0 * nHepp
  
  dnH0_dt = k2(T) * nHp * ne - k1(T) * nH0 * ne
  
  dnHe0_dt = k5(T) * nHep * ne + k7(T) * nHep * ne - k3(T) * nHe0 * ne
  
  dnHep_dt = k6(T) * nHepp * ne + k3(T) * nHe0 * ne - k4(T) * nHep * ne - k5(T) * nHep * ne - k7(T) * nHep * ne
  
  Lamb = Lambda(T, nH, nH0, nHe0, nHep)
  
  ntot = n_tot(nH, nH0, nHe0, nHep)

  dT_dt = -1.0 * (gamma - 1.0) / kB / ntot * Lamb
  
  return np.array([dnH0_dt, dnHe0_dt, dnHep_dt, dT_dt])




#----- p_interp
def g_interp(Tref, T0, gx):

  nt = -1
  
  N = len(Tref)
  
  for j in range(N):
    
    if (Tref[j] >= T0) & (nt == -1):
      nt = j
  
  if nt == -1:
    nt = N - 1 
  
  diff = Tref[nt] - Tref[nt - 1]
  fhi = (T0 - Tref[nt - 1]) / diff
  flow = 1.0 - fhi
  
  gnew = flow * np.log10(gx[nt - 1]) + fhi * np.log10(gx[nt])
  
  return 10**gnew
  




gamma = 5./3.
kB = 1.3807e-16
X = 0.76
Y = 1.0 - X

nH = 100.0

with open('HeH_cooling_rates.pkl', 'rb') as f:
  dictx = pickle.load(f)

Tref = dictx['Tref']

g1 = dictx['g1']
g2 = dictx['g2']
g3 = dictx['g3']
g4 = dictx['g4']
g5 = dictx['g5']
g6 = dictx['g6']
g7 = dictx['g7']
g8 = dictx['g8']
g9 = dictx['g9']
g10= dictx['g10']


He_solar = 10**(-1.07)
nHe = He_solar * nH
print('nHe (cm^-3) = ', nHe)


y0 = np.array([1e-3, 1.6e-6, 6e-3, 1e6])

num_time_pts = 3000
dt = 0.01 * 3.16e7

print(f'System will evolve for {(num_time_pts*dt/3.16e7):.2f} years!')

t = 0.0

Y = np.zeros((4, num_time_pts))
Y[:, 0] = y0
yin = y0

t_yrs = np.zeros(num_time_pts)

for i in range(num_time_pts - 1):

  yout = rk4SingleStep(func, dt, t, yin)
  Y[:, i+1] = yout
  yin = yout
  
  t_yrs[i+1] = t
  
  t += dt


print(Y)

T = Y[-1, :]

nH0 = Y[0, :]
nHp = nH - nH0
nHe0 = Y[1, :]
nHep = Y[2, :]
nHepp = nHe - nHe0 - nHep


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

plt.figure(figsize = (16, 14))

plt.subplot(2, 3, 1)
plt.scatter(t_yrs/3.16e7, np.log10(T), s = 5, color = 'k', label = 'my own code')
plt.scatter(t_Arr_in_yrsx, np.log10(TEvolx), s = 2, color = 'orange', label = 'chimes result', linestyle = '--')
plt.xlim(0, 3000)
plt.ylim(3, 8)
plt.legend()

plt.subplot(2, 3, 2)

plt.plot(t_yrs/3.16e7, nHe0/nHe, color = 'r', label = 'nHe0')
plt.plot(t_yrs/3.16e7, nHep/nHe, color = 'g', label = 'nHep')
plt.plot(t_yrs/3.16e7, nHepp/nHe, color = 'b', label = 'nHepp')

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

plt.savefig('rk4_HeH.png')

plt.show()





