import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import interp1d


#-----------------------------------------------------------
#------------> Cooling rates from CHIMES table <------------
#-----------------------------------------------------------
with h5py.File('chimes_main_data.hdf5', 'r') as file:
  
  Temp_4d = file['TableBins/cool_4d_Temperatures'][:]
  HIDensity_4d = file['TableBins/cool_4d_HIDensities'][:] # only used for low T
  elecDensity_4d = file['TableBins/cool_4d_ElectronDensities'][:] # only used for low T
  HIIDensity_4d = file['TableBins/cool_4d_HIIDensities'][:] # only used for low T
  
  rates_4d = file['cooling/rates_4d'][:] # NOTE it is rates_4d for low Temp
  
  #-------- hiT_4d ---------------
  Temp_hiT_4d = file['TableBins/cool_hiT_4d_Temperatures'][:]
  rates_hiT_4d = file['cooling/rates_hiT_4d'][:] # NOTE it is rates_4d for high Temp
  
  #------------------------------------------------------------
  #------ CII section -----------------------------------------
  #------------------------------------------------------------
  coolants_2d = file['cooling/coolants_2d'][:] # is the same for low and hiT.... Only Temp will change!! And in hiT it is only a function of T!!!
  
  Temp_2d = file['TableBins/cool_2d_Temperatures'][:]
  elecDensity_2d = file['TableBins/cool_2d_ElectronDensities'][:]
  
  rates_2d = file['cooling/rates_2d'][:] # NOTE it is rates_2d for low Temp
  
  #-------- hiT_2d ---------------
  Temp_hiT_2d = file['TableBins/cool_hiT_2d_Temperatures'][:]
  rates_hiT_2d = file['cooling/rates_hiT_2d'][:] # NOTE it is rates_2d for high Temp


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
  else:
    C0_rates = rates_hiT_4d[0, :]
    interp_4d = interp1d(Temp_hiT_4d, C0_rates, kind='linear', fill_value="extrapolate")
    res = interp_4d(T)

  return res


#----- Cp_cooling_rate 
def Cp_cooling_rate(T, nelec, Temp_2d, elecDensity_2d):

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


R_low = C0_cooling_rate(10**3.0, 10**-8.0, 10**-8.0, 10**-8.0, Temp_4d, HIDensity_4d, elecDensity_4d, HIIDensity_4d)
print('C0: Rlow = ', R_low)

R_hi = C0_cooling_rate(10**5.15, 10**-8.0, 10**-8.0, 10**-8.0, Temp_4d, HIDensity_4d, elecDensity_4d, HIIDensity_4d)
print('C0: Rhi = ', R_hi)

test_lowT = Cp_cooling_rate(10**3.0, 10**-8.0, Temp_2d, elecDensity_2d)
print('C1: lowT test = ', test_lowT)


test_hiT_T = Cp_cooling_rate(10**5.15, 10**-8.0, Temp_2d, elecDensity_2d)
print('C1: hiT test = ', test_hiT_T)






