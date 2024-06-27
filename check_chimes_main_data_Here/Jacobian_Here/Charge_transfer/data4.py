
# In this version, I include grain_recombination (29 May 2024).

import h5py
from scipy.interpolate import interp1d
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import interp2d

gamma = 5./3.
kB = 1.3807e-16

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
  #-------- C0_cooling_ratehiT_2d ---------------
  Temp_hiT_2d = file['TableBins/cool_hiT_2d_Temperatures'][:]
  rates_hiT_2d = file['cooling/rates_hiT_2d'][:] # NOTE it is rates_2d for high Temp

  #----------> grain_recombination <-------------
  grain_recomb_rates = file['grain_recombination/rates'][:]
  Psi = file['TableBins/Psi'][:]
  grain_cooling_rates = file['cooling/grain_recombination'][:]
  grain_cooling_rates = 10**grain_cooling_rates

# NOTE: I used "T_dependent_reactants.py" code to find indices 111, 0, 108, etc!!!
#---- Reaction rates -------
k1x = rates[111, :] # Reaction: (H0 + e ---> Hp + 2e) ::: H0 Collisional ionization
k2x = ratesAB[0, :] # Reaction: (Hp + e ---> H0 + γ)  ::: photo-recombination of Hp and e.
k3x = rates[108, :] # Reaction: (He0 + e ---> Hep + 2e) ::: He0 Collisional ionization 
k4x = rates[0, :]   # Reaction: (Hep + e ---> Hepp + 2e) ::: Hep Collisional ionization
k5x = ratesAB[1, :] # Reaction: (Hep + e ---> He0 + γ)  ::: photo-recombination of Hep and e. !!! PROBABLY di-electric is also included ?????
k6x = rates[222, :] # Reaction: (Hepp + e ---> Hep + γ)  ::: photo-recombination of Hepp and e.

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

ct_C_0_1_Hepx = rates[218, :] # C0 + Hep ---> C1 + He0
ct_C_0_1_Hpx  = rates[219, :] # C0 + Hp  ---> C1 + H0

ct_C_1_0_H0x  = rates[214, :] # C1 + H0  ---> C0 + Hp
ct_C_1_2_Hepx = rates[215, :] # C1 + Hep ---> C2 + He0 

ct_C_2_1_H0x  = rates[211, :] # C2 + H0 ---> C1 + Hp

ct_C_3_2_He0x = rates[207, :] # C3 + He0 ---> C2 + Hep
ct_C_3_2_H0x  = rates[208, :] # C3 + H0  ---> C2 + Hp

ct_C_4_3_H0x  = rates[225, :] # C4 + H0  ---> C3 + Hp
ct_C_4_3_He0x = rates[244, :] # C4 + He0 ---> C3 + Hep

ct_C_5_4_H0x  = rates[242, :] # C5 + H0  ---> C4 + Hp

#---- Cooling rates ------
g1x = cooling_rates[0, :] # cooling via H0
g2x = cooling_rates[1, :] # cooling via Hp
g3x = cooling_rates[2, :] # cooling via He0
g4x = cooling_rates[3, :] # cooling via Hep
g5x = cooling_rates[4, :] # cooling via Hepp

# cooling via C0 and C1 are determined using "C0_cooling_rate" and "C0_cooling_rate" functions inside the "Lambda" function directly!!
gC2x = cooling_rates[5, :]
gC3x = cooling_rates[6, :]
gC4x = cooling_rates[7, :]
gC5x = cooling_rates[8, :]
gC6x = cooling_rates[9, :]


k1 = interp1d(Temp, k1x, kind='linear', fill_value="extrapolate")
k2 = interp1d(Temp, k2x, kind='linear', fill_value="extrapolate")
k3 = interp1d(Temp, k3x, kind='linear', fill_value="extrapolate")
k4 = interp1d(Temp, k4x, kind='linear', fill_value="extrapolate")
k5 = interp1d(Temp, k5x, kind='linear', fill_value="extrapolate")
k6 = interp1d(Temp, k6x, kind='linear', fill_value="extrapolate")

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


ct_C_0_1_Hep = interp1d(Temp, ct_C_0_1_Hepx, kind='linear', fill_value="extrapolate") # C0 + Hep ---> C1 + He0
ct_C_0_1_Hp  = interp1d(Temp, ct_C_0_1_Hpx, kind='linear', fill_value="extrapolate") # C0 + Hp  ---> C1 + H0

ct_C_1_0_H0  = interp1d(Temp, ct_C_1_0_H0x, kind='linear', fill_value="extrapolate") # C1 + H0  ---> C0 + Hp
ct_C_1_2_Hep = interp1d(Temp, ct_C_1_2_Hepx, kind='linear', fill_value="extrapolate") # C1 + Hep ---> C2 + He0 

ct_C_2_1_H0  = interp1d(Temp, ct_C_2_1_H0x, kind='linear', fill_value="extrapolate") # C2 + H0 ---> C1 + Hp

ct_C_3_2_He0 = interp1d(Temp, ct_C_3_2_He0x, kind='linear', fill_value="extrapolate") # C3 + He0 ---> C2 + Hep
ct_C_3_2_H0  = interp1d(Temp, ct_C_3_2_H0x, kind='linear', fill_value="extrapolate") # C3 + H0  ---> C2 + Hp

ct_C_4_3_H0  = interp1d(Temp, ct_C_4_3_H0x, kind='linear', fill_value="extrapolate") # C4 + H0  ---> C3 + Hp
ct_C_4_3_He0 = interp1d(Temp, ct_C_4_3_He0x, kind='linear', fill_value="extrapolate") # C4 + He0 ---> C3 + Hep

ct_C_5_4_H0  = interp1d(Temp, ct_C_5_4_H0x, kind='linear', fill_value="extrapolate") # C5 + H0  ---> C4 + Hp


g1 = interp1d(Temp, g1x, kind='linear', fill_value="extrapolate")
g2 = interp1d(Temp, g2x, kind='linear', fill_value="extrapolate")
g3 = interp1d(Temp, g3x, kind='linear', fill_value="extrapolate")
g4 = interp1d(Temp, g4x, kind='linear', fill_value="extrapolate")
g5 = interp1d(Temp, g5x, kind='linear', fill_value="extrapolate")

gC2 = interp1d(Temp, gC2x, kind='linear', fill_value="extrapolate")
gC3 = interp1d(Temp, gC3x, kind='linear', fill_value="extrapolate")
gC4 = interp1d(Temp, gC4x, kind='linear', fill_value="extrapolate")
gC5 = interp1d(Temp, gC5x, kind='linear', fill_value="extrapolate")
gC6 = interp1d(Temp, gC6x, kind='linear', fill_value="extrapolate")


#----- C0_cooling_rate 
def C0_cooling_rate(T, nHI, nelec, nHII, Temp_4d, HIDensity_4d, elecDensity_4d, HIIDensity_4d):

  T = np.log10(T)
  nHI = np.log10(nHI)
  nelec = np.log10(nelec)
  nHII = np.log10(nHII)

  if T <= 4:
    C0_rates = rates_4d[0, :]
    interp_4d = RegularGridInterpolator((Temp_4d, HIDensity_4d, elecDensity_4d, HIIDensity_4d), C0_rates)
    res = interp_4d(np.array([T, nHI, nelec, nHII]))[0]
    #print(T, nelec, np.log10(10**res / 10**nelec))
    res = np.log10(10**res / 10**nelec) # because ne is multiplied later so should be cancelled by this. Intrinsically ne was included in chimes Table!!!
  else:
    C0_rates = rates_hiT_4d[0, :]
    interp_4d = interp1d(Temp_hiT_4d, C0_rates, kind='linear', fill_value="extrapolate")
    res = interp_4d(T)
    #print(T, nelec, res)

  return res


#----- Cp_cooling_rate 
def Cp_cooling_rate(T, nelec, Temp_2d, elecDensity_2d): # include Temp_hiT here !!!!!!!!!!!!!!!!!!!!!!!

  T = np.log10(T)
  nelec = np.log10(nelec)

  if T <= 4.95:
    Cp_rates = rates_2d[0, :]
    interp_2d = RegularGridInterpolator((Temp_2d, elecDensity_2d), Cp_rates)
    res = interp_2d(np.array([T, nelec]))[0]
    #print(T, nelec, res, np.log10(10**res / 10**nelec))
    #res = np.log10(10**res / 10**nelec)
  else:
    Cp_rates = rates_hiT_2d[0, :]
    interp_2d = interp1d(Temp_hiT_2d, Cp_rates, kind='linear', fill_value="extrapolate")
    res = interp_2d(T)
    #print(T, nelec, res)

  return res



#----- grain_recomb_rate
def grain_recomb_rate(ionx, T, nelec, G0, A_v, Temp, Psi): # Note: Temp and Psi are in log from CHIMES table!
  
  if ionx == 'Hp':
    i = 0
  elif ionx == 'Hep':
    i = 1
  elif ionx == 'C1':
    i = 2
  
  Psix = G0 * np.exp(-2.77 * A_v) * T**0.5 / nelec
  Psix = np.log10(Psix)
  
  T = np.log10(T)
  
  interp_2d = RegularGridInterpolator((Temp, Psi), grain_recomb_rates[i, :])
  res = interp_2d(np.array([T, Psix]))[0]

  return res



#----- grain_cool_rate
def grain_cool_rate(T, nelec, nH0, nHp, G0, A_v, dust_ratio, Temp, Psi): # Note: Temp and Psi are in log from CHIMES table!
  
  Psix = G0 * np.exp(-2.77 * A_v) * T**0.5 / nelec
  Psix = np.log10(Psix)
  
  nHtot = nH0 + nHp
  
  cool_rates = grain_cooling_rates * dust_ratio * nelec / nHtot # Note: none of grain_cooling_rates, nelec and nHtot is in log form!!
  cool_rates = np.log10(cool_rates)
  
  interp_2d = RegularGridInterpolator((Temp, Psi), cool_rates)
  
  T = np.log10(T)
  
  res = interp_2d(np.array([T, Psix]))[0]

  #print(T, res)

  return res









