
import h5py
from scipy.interpolate import interp1d
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import interp2d
import time


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



R_CI_to_CII_via_HeII_ = rates[218, :]
R_CI_to_CII_via_HeII = interp1d(Temp, R_CI_to_CII_via_HeII_, kind="linear", fill_value="extrapolate")


print(rates.shape)
print(Temp.shape)

#print(Temp)


T = T0 = np.log10(20000.)


#----- find_index
def find_index(Temp, T0):

  N_T = len(Temp)

  for i in range(N_T):
    if Temp[i] >= T0:
      ndx_T = i
      break

  ndxL = i

  if ndx_T == N_T - 1:
    ndxL = N_T - 2

  if ndx_T == 0:
    ndxL = 0
  
  return ndxL # ndxR will be created outside by adding 1 to ndxL !





print(R_CI_to_CII_via_HeII_[ndxL])
print(R_CI_to_CII_via_HeII_[ndxR])


#
T1 = time.time()
RR = R_CI_to_CII_via_HeII(T)
print(f'Elapsed time = {time.time() - T1}')
print(RR)






