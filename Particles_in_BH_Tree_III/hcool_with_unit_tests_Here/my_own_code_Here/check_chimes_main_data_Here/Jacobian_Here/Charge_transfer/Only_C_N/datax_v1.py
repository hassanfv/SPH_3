
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

k1 = interp1d(Temp, k1x, kind='linear', fill_value="extrapolate")
k2 = interp1d(Temp, k2x, kind='linear', fill_value="extrapolate")
k3 = interp1d(Temp, k3x, kind='linear', fill_value="extrapolate")
k4 = interp1d(Temp, k4x, kind='linear', fill_value="extrapolate")
k5 = interp1d(Temp, k5x, kind='linear', fill_value="extrapolate")
k6 = interp1d(Temp, k6x, kind='linear', fill_value="extrapolate")

R_CI_to_CII_via_HeII_ = rates[218, :]
R_CI_to_CII_via_HII_ = rates[219, :]
R_CI_to_CII_via_e_ = rates[220, :]
R_CII_to_CI_via_HI_ = rates[214, :]
R_CII_to_CIII_via_HeII_ = rates[215, :]
R_CII_to_CI_via_e_ = rates[216, :]
R_CII_to_CIII_via_e_ = rates[217, :]
R_CIII_to_CII_via_HI_ = rates[211, :]
R_CIII_to_CII_via_e_ = rates[212, :]
R_CIII_to_CIV_via_e_ = rates[213, :]
R_CIV_to_CIII_via_HeI_ = rates[207, :]
R_CIV_to_CIII_via_HI_ = rates[208, :]
R_CIV_to_CIII_via_e_ = rates[209, :]
R_CIV_to_CV_via_e_ = rates[210, :]
R_CV_to_CVI_via_e_ = rates[206, :]
R_CV_to_CIV_via_e_ = rates[223, :]
R_CV_to_CIV_via_HI_ = rates[225, :]
R_CV_to_CIV_via_HeI_ = rates[244, :]
R_CVI_to_CVII_via_e_ = rates[226, :]
R_CVI_to_CV_via_HI_ = rates[242, :]
R_CVI_to_CV_via_e_ = rates[243, :]
R_CVII_to_CVI_via_e_ = rates[241, :]
R_Cm_to_CI_via_HII_ = rates[105, :]

R_CI_to_CII_via_HeII = interp1d(Temp, R_CI_to_CII_via_HeII_, kind="linear", fill_value="extrapolate")
R_CI_to_CII_via_HII = interp1d(Temp, R_CI_to_CII_via_HII_, kind="linear", fill_value="extrapolate")
R_CI_to_CII_via_e = interp1d(Temp, R_CI_to_CII_via_e_, kind="linear", fill_value="extrapolate")
R_CII_to_CI_via_HI = interp1d(Temp, R_CII_to_CI_via_HI_, kind="linear", fill_value="extrapolate")
R_CII_to_CIII_via_HeII = interp1d(Temp, R_CII_to_CIII_via_HeII_, kind="linear", fill_value="extrapolate")
R_CII_to_CI_via_e = interp1d(Temp, R_CII_to_CI_via_e_, kind="linear", fill_value="extrapolate")
R_CII_to_CIII_via_e = interp1d(Temp, R_CII_to_CIII_via_e_, kind="linear", fill_value="extrapolate")
R_CIII_to_CII_via_HI = interp1d(Temp, R_CIII_to_CII_via_HI_, kind="linear", fill_value="extrapolate")
R_CIII_to_CII_via_e = interp1d(Temp, R_CIII_to_CII_via_e_, kind="linear", fill_value="extrapolate")
R_CIII_to_CIV_via_e = interp1d(Temp, R_CIII_to_CIV_via_e_, kind="linear", fill_value="extrapolate")
R_CIV_to_CIII_via_HeI = interp1d(Temp, R_CIV_to_CIII_via_HeI_, kind="linear", fill_value="extrapolate")
R_CIV_to_CIII_via_HI = interp1d(Temp, R_CIV_to_CIII_via_HI_, kind="linear", fill_value="extrapolate")
R_CIV_to_CIII_via_e = interp1d(Temp, R_CIV_to_CIII_via_e_, kind="linear", fill_value="extrapolate")
R_CIV_to_CV_via_e = interp1d(Temp, R_CIV_to_CV_via_e_, kind="linear", fill_value="extrapolate")
R_CV_to_CVI_via_e = interp1d(Temp, R_CV_to_CVI_via_e_, kind="linear", fill_value="extrapolate")
R_CV_to_CIV_via_e = interp1d(Temp, R_CV_to_CIV_via_e_, kind="linear", fill_value="extrapolate")
R_CV_to_CIV_via_HI = interp1d(Temp, R_CV_to_CIV_via_HI_, kind="linear", fill_value="extrapolate")
R_CV_to_CIV_via_HeI = interp1d(Temp, R_CV_to_CIV_via_HeI_, kind="linear", fill_value="extrapolate")
R_CVI_to_CVII_via_e = interp1d(Temp, R_CVI_to_CVII_via_e_, kind="linear", fill_value="extrapolate")
R_CVI_to_CV_via_HI = interp1d(Temp, R_CVI_to_CV_via_HI_, kind="linear", fill_value="extrapolate")
R_CVI_to_CV_via_e = interp1d(Temp, R_CVI_to_CV_via_e_, kind="linear", fill_value="extrapolate")
R_CVII_to_CVI_via_e = interp1d(Temp, R_CVII_to_CVI_via_e_, kind="linear", fill_value="extrapolate")
R_Cm_to_CI_via_HII = interp1d(Temp, R_Cm_to_CI_via_HII_, kind="linear", fill_value="extrapolate")







