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



R_low = C0_cooling_rate(10**3.0, 10**-8.0, 10**-8.0, 10**-8.0, Temp_4d, HIDensity_4d, elecDensity_4d, HIIDensity_4d)
print('Rlow = ', R_low)

R_hi = C0_cooling_rate(10**5.15, 10**-8.0, 10**-8.0, 10**-8.0, Temp_4d, HIDensity_4d, elecDensity_4d, HIIDensity_4d)
print('Rhi = ', R_hi)


