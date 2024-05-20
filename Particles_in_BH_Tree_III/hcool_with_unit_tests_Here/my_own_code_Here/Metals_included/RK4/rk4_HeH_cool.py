
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





#----------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
#------------------------------------- Cooling Rates (erg.s^-1.cm^3) --------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
# Cooling rate due to H collisional ionization: (H0 + e ---> Hp + 2e)  ::: in erg.s^-1.cm^3  ::: will be multiplied by nH0 and ne later in the code.
def g1(T):

  g1_val = 1.27e-21 * T**0.5 * np.exp(-157809.1/T) / (1.0 + (T/1e5)**0.5)
  
  return g1_val


# Cooling rate due to Hp photo-recombination: (Hp + e ---> H0 + γ)  ::: in erg.s^-1.cm^3  ::: will be multiplied by nHp and ne later in the code.
def g2(T):

  #g2_val = 8.70e-27 * T**0.5 / (T/1e3)**0.2 / (1.0 + (T/1e6)**0.7)
  g2_val = kB * T * k2(T) # See LaMothe and Ferland - 2001, page 1, the lines below equation 1!

  return g2_val


# Cooling due to H0 collisional excitation via collision with electrons ::: in erg.s^-1.cm^3 ::: will be multiplied by nH0 and ne later in the code.
def g3(T):

  g3_val = 7.50e-19 * np.exp(-118348/T) / (1.0 + (T/1e5)**0.5)
  
  return g3_val


# Cooling due to free-free emission ::: will be multiplied by (nHp+nHep+nHepp)*ne later in the code.
def g4(T):

  gff = 1.1 + 0.34 * np.exp(-(5.5 - np.log10(T))**2/3.0)
  g4_val = 1.42e-27 * gff * T**0.5
  
  return g4_val


# Cooling due to Hep collisional excitation ::: will be multiplied by nHep and ne later in the code
def g5(T):

  g5_val = 5.54e-17 / T**0.397 * np.exp(-473638./T) / (1. + (T/1e5)**0.5)
  
  return g5_val


# Cooling due to He0 collisional ioniozation ::: will be multiplied by nHe0 and ne later in the code
def g6(T):
  
  g6_val = 9.38e-22 * T**0.5 * np.exp(-285335.4/T) / (1. + (T/1e5)**0.5)
  
  return g6_val


# Cooling via Hep collisional ionization ::: will be multiplied by nHep and ne later in the code
def g7(T):

  g7_val = 4.95e-22 * T**0.5 * np.exp(-631515./T) / (1. + (T/1e5)**0.5)
  
  return g7_val


# Cooling via Hep recombination with electron ::: will be multiplied by nHep and ne later in the code
def g8(T):

  g8_val = 1.55e-26 * T**0.3647
  
  return g8_val


# Cooling via Hepp recombination with electron ::: will be multiplied by nHepp and ne later in the code
def g9(T):
  
  g9_val = 3.48e-26 * T**0.5 / (T/1e3)**0.2 / (1. + (T/1e6)**0.7)
  
  return g9_val


# Cooling via di-electronic recombination of Hep with electrons ::: will be multiplied by nHep and ne later in the code
def g10(T):
  
  g10_val = 1.24e-13 / T**1.5 * np.exp(-470000./T) * (1. + 0.3 * np.exp(-94000./T))

  return g10_val





# Total cooling rate in erg.s^-1.cm^-3 ===>NOTE that Hep and Hepp is excluded in free-free here as we are only considering H in this code!
def Lambda(T, nH, nH0, nHe0, nHep):

  nHp = nH - nH0
  nHepp = nHe - nHe0 - nHep
  
  ne = nHp + nHep + 2.0 * nHepp
  
  Lamb = (
           g1(T) * ne * nH0 # collisional ionization of hydrogen ::::::: (H0 + e ---> Hp + 2e).
         + g2(T)  * ne * nHp # photo-recombination of Hp with electron :: (Hp + e ---> H0 + γ).
         + g3(T)  * ne * nH0 # collisional excitaion of H0.
         + g4(T) * ne * (nHp + nHep + 4.0 * nHepp) # free-free emission
         + g5(T) * nHep * ne # collisional excitation of Hep.
         + g6(T) * nHe0 * ne # He0 collisional ionization
         + g7(T) * nHep * ne # Hep collisional ionization
         + g8(T) * nHep * ne # Hep recombination to He0
         + g9(T) * nHepp * ne# Hepp recombination to Hep
         + g10(T) * nHep * ne# Hep di-electric recombination to He0
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



gamma = 5./3.
kB = 1.3807e-16
X = 0.76
Y = 1.0 - X

nH = 1000.0

He_solar = 10**(-1.07)
nHe = He_solar * nH
print('nHe (cm^-3) = ', nHe)


y0 = np.array([1e-3, 1.6e-6, 6e-3, 1e6])

num_time_pts = 600000
dt = 0.005 * 3.16e7

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





