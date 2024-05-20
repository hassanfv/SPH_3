import h5py
from scipy.interpolate import interp1d
import numpy as np
from scipy.interpolate import RegularGridInterpolator

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
    #return -30
    C0_rates = rates_4d[0, :]
    interp_4d = RegularGridInterpolator((Temp_4d, HIDensity_4d, elecDensity_4d, HIIDensity_4d), C0_rates)
    res = interp_4d(np.array([T, nHI, nelec, nHII]))[0]
    res = np.log10(10**res / 10**nelec) # because ne is added later so should be cancelled by this. Intrinsically ne included!!!
  else:
    C0_rates = rates_hiT_4d[0, :]
    interp_4d = interp1d(Temp_hiT_4d, C0_rates, kind='linear', fill_value="extrapolate")
    res = interp_4d(T)

  return res


#----- Cp_cooling_rate 
def Cp_cooling_rate(T, nelec, Temp_2d, elecDensity_2d): # include Temp_hiT here !!!!!!!!!!!!!!!!!!!!!!!

  T = np.log10(T)
  nelec = np.log10(nelec)

  if T <= 4:
    return -30
    #Cp_rates = rates_2d[0, :]
    #interp_2d = RegularGridInterpolator((Temp_2d, elecDensity_2d), Cp_rates)
    #res = interp_2d(np.array([T, nelec]))[0]
  else:
    Cp_rates = rates_hiT_2d[0, :]
    interp_2d = interp1d(Temp_hiT_2d, Cp_rates, kind='linear', fill_value="extrapolate")
    res = interp_2d(T)

  return res





