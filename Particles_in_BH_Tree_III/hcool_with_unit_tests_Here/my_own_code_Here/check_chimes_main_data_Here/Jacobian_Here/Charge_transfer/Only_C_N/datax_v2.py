
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

R_HII_to_HI_via_e_caseA_ = ratesAB[0, :] # H CaseA
R_HeII_to_HeI_via_e_caseA_ = ratesAB[1, :] # He CaseA
R_HI_to_HII_via_e_ = rates[111, :]
R_HI_to_Hm_via_e_ = rates[113, :]
R_Hm_to_HI_via_HI_ = rates[109, :]
R_Hm_to_HI_via_e_ = rates[110, :]
R_Hm_to_HI_via_HII_ = rates[112, :]
R_HeI_to_HeII_via_HII_ = rates[106, :]
R_HeI_to_HeII_via_e_ = rates[108, :]
R_HeII_to_HeIII_via_e_ = rates[0, :]
R_HeII_to_HeI_via_Hm_ = rates[103, :]
R_HeII_to_HeI_via_HI_ = rates[107, :]
R_HeIII_to_HeII_via_HI_ = rates[221, :]
R_HeIII_to_HeII_via_e_ = rates[222, :]

R_HII_to_HI_via_e_caseA = interp1d(Temp, R_HII_to_HI_via_e_caseA_, kind="linear", fill_value="extrapolate") # H CaseA
R_HeII_to_HeI_via_e_caseA = interp1d(Temp, R_HeII_to_HeI_via_e_caseA_, kind="linear", fill_value="extrapolate") # He CaseA
R_HI_to_HII_via_e = interp1d(Temp, R_HI_to_HII_via_e_, kind="linear", fill_value="extrapolate")
R_HI_to_Hm_via_e = interp1d(Temp, R_HI_to_Hm_via_e_, kind="linear", fill_value="extrapolate")
R_Hm_to_HI_via_HI = interp1d(Temp, R_Hm_to_HI_via_HI_, kind="linear", fill_value="extrapolate")
R_Hm_to_HI_via_e = interp1d(Temp, R_Hm_to_HI_via_e_, kind="linear", fill_value="extrapolate")
R_Hm_to_HI_via_HII = interp1d(Temp, R_Hm_to_HI_via_HII_, kind="linear", fill_value="extrapolate")
R_HeI_to_HeII_via_HII = interp1d(Temp, R_HeI_to_HeII_via_HII_, kind="linear", fill_value="extrapolate")
R_HeI_to_HeII_via_e = interp1d(Temp, R_HeI_to_HeII_via_e_, kind="linear", fill_value="extrapolate")
R_HeII_to_HeIII_via_e = interp1d(Temp, R_HeII_to_HeIII_via_e_, kind="linear", fill_value="extrapolate")
R_HeII_to_HeI_via_Hm = interp1d(Temp, R_HeII_to_HeI_via_Hm_, kind="linear", fill_value="extrapolate")
R_HeII_to_HeI_via_HI = interp1d(Temp, R_HeII_to_HeI_via_HI_, kind="linear", fill_value="extrapolate")
R_HeIII_to_HeII_via_HI = interp1d(Temp, R_HeIII_to_HeII_via_HI_, kind="linear", fill_value="extrapolate")
R_HeIII_to_HeII_via_e = interp1d(Temp, R_HeIII_to_HeII_via_e_, kind="linear", fill_value="extrapolate")


#---- Cooling rates ------
g1x = cooling_rates[0, :] # cooling via H0
g2x = cooling_rates[1, :] # cooling via Hp
g3x = cooling_rates[2, :] # cooling via He0
g4x = cooling_rates[3, :] # cooling via Hep
g5x = cooling_rates[4, :] # cooling via Hepp

g1 = interp1d(Temp, g1x, kind='linear', fill_value="extrapolate")
g2 = interp1d(Temp, g2x, kind='linear', fill_value="extrapolate")
g3 = interp1d(Temp, g3x, kind='linear', fill_value="extrapolate")
g4 = interp1d(Temp, g4x, kind='linear', fill_value="extrapolate")
g5 = interp1d(Temp, g5x, kind='linear', fill_value="extrapolate")


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


gCIII_ = cooling_rates[5, :]
gCIV_ = cooling_rates[6, :]
gCV_ = cooling_rates[7, :]
gCVI_ = cooling_rates[8, :]
gCVII_ = cooling_rates[9, :]

gCIII = interp1d(Temp, gCIII_, kind="linear", fill_value="extrapolate")
gCIV = interp1d(Temp, gCIV_, kind="linear", fill_value="extrapolate")
gCV = interp1d(Temp, gCV_, kind="linear", fill_value="extrapolate")
gCVI = interp1d(Temp, gCVI_, kind="linear", fill_value="extrapolate")
gCVII = interp1d(Temp, gCVII_, kind="linear", fill_value="extrapolate")


#----------- Nitrogen Section --------------
#--- REACTION RATES
R_NI_to_NII_via_HII_ = rates[239, :]
R_NI_to_NII_via_e_ = rates[240, :]
R_NII_to_NI_via_HI_ = rates[235, :]
R_NII_to_NIII_via_HeII_ = rates[236, :]
R_NII_to_NI_via_e_ = rates[237, :]
R_NII_to_NIII_via_e_ = rates[238, :]
R_NIII_to_NII_via_HeI_ = rates[231, :]
R_NIII_to_NII_via_HI_ = rates[232, :]
R_NIII_to_NII_via_e_ = rates[233, :]
R_NIII_to_NIV_via_e_ = rates[234, :]
R_NIV_to_NIII_via_HeI_ = rates[227, :]
R_NIV_to_NIII_via_HI_ = rates[228, :]
R_NIV_to_NIII_via_e_ = rates[229, :]
R_NIV_to_NV_via_e_ = rates[230, :]
R_NV_to_NIV_via_HeI_ = rates[202, :]
R_NV_to_NIV_via_HI_ = rates[203, :]
R_NV_to_NIV_via_e_ = rates[204, :]
R_NV_to_NVI_via_e_ = rates[205, :]
R_NVI_to_NV_via_HI_ = rates[178, :]
R_NVI_to_NVII_via_e_ = rates[180, :]
R_NVI_to_NV_via_e_ = rates[246, :]
R_NVII_to_NVI_via_e_ = rates[176, :]
R_NVII_to_NVIII_via_e_ = rates[177, :]
R_NVIII_to_NVII_via_e_ = rates[175, :]

R_NI_to_NII_via_HII = interp1d(Temp, R_NI_to_NII_via_HII_, kind="linear", fill_value="extrapolate")
R_NI_to_NII_via_e = interp1d(Temp, R_NI_to_NII_via_e_, kind="linear", fill_value="extrapolate")
R_NII_to_NI_via_HI = interp1d(Temp, R_NII_to_NI_via_HI_, kind="linear", fill_value="extrapolate")
R_NII_to_NIII_via_HeII = interp1d(Temp, R_NII_to_NIII_via_HeII_, kind="linear", fill_value="extrapolate")
R_NII_to_NI_via_e = interp1d(Temp, R_NII_to_NI_via_e_, kind="linear", fill_value="extrapolate")
R_NII_to_NIII_via_e = interp1d(Temp, R_NII_to_NIII_via_e_, kind="linear", fill_value="extrapolate")
R_NIII_to_NII_via_HeI = interp1d(Temp, R_NIII_to_NII_via_HeI_, kind="linear", fill_value="extrapolate")
R_NIII_to_NII_via_HI = interp1d(Temp, R_NIII_to_NII_via_HI_, kind="linear", fill_value="extrapolate")
R_NIII_to_NII_via_e = interp1d(Temp, R_NIII_to_NII_via_e_, kind="linear", fill_value="extrapolate")
R_NIII_to_NIV_via_e = interp1d(Temp, R_NIII_to_NIV_via_e_, kind="linear", fill_value="extrapolate")
R_NIV_to_NIII_via_HeI = interp1d(Temp, R_NIV_to_NIII_via_HeI_, kind="linear", fill_value="extrapolate")
R_NIV_to_NIII_via_HI = interp1d(Temp, R_NIV_to_NIII_via_HI_, kind="linear", fill_value="extrapolate")
R_NIV_to_NIII_via_e = interp1d(Temp, R_NIV_to_NIII_via_e_, kind="linear", fill_value="extrapolate")
R_NIV_to_NV_via_e = interp1d(Temp, R_NIV_to_NV_via_e_, kind="linear", fill_value="extrapolate")
R_NV_to_NIV_via_HeI = interp1d(Temp, R_NV_to_NIV_via_HeI_, kind="linear", fill_value="extrapolate")
R_NV_to_NIV_via_HI = interp1d(Temp, R_NV_to_NIV_via_HI_, kind="linear", fill_value="extrapolate")
R_NV_to_NIV_via_e = interp1d(Temp, R_NV_to_NIV_via_e_, kind="linear", fill_value="extrapolate")
R_NV_to_NVI_via_e = interp1d(Temp, R_NV_to_NVI_via_e_, kind="linear", fill_value="extrapolate")
R_NVI_to_NV_via_HI = interp1d(Temp, R_NVI_to_NV_via_HI_, kind="linear", fill_value="extrapolate")
R_NVI_to_NVII_via_e = interp1d(Temp, R_NVI_to_NVII_via_e_, kind="linear", fill_value="extrapolate")
R_NVI_to_NV_via_e = interp1d(Temp, R_NVI_to_NV_via_e_, kind="linear", fill_value="extrapolate")
R_NVII_to_NVI_via_e = interp1d(Temp, R_NVII_to_NVI_via_e_, kind="linear", fill_value="extrapolate")
R_NVII_to_NVIII_via_e = interp1d(Temp, R_NVII_to_NVIII_via_e_, kind="linear", fill_value="extrapolate")
R_NVIII_to_NVII_via_e = interp1d(Temp, R_NVIII_to_NVII_via_e_, kind="linear", fill_value="extrapolate")

#--- COOLING RATES
gNI_ = cooling_rates[10, :]
gNIII_ = cooling_rates[11, :]
gNIV_ = cooling_rates[12, :]
gNV_ = cooling_rates[13, :]
gNVI_ = cooling_rates[14, :]
gNVII_ = cooling_rates[15, :]
gNVIII_ = cooling_rates[16, :]

gNI = interp1d(Temp, gNI_, kind="linear", fill_value="extrapolate")
gNIII = interp1d(Temp, gNIII_, kind="linear", fill_value="extrapolate")
gNIV = interp1d(Temp, gNIV_, kind="linear", fill_value="extrapolate")
gNV = interp1d(Temp, gNV_, kind="linear", fill_value="extrapolate")
gNVI = interp1d(Temp, gNVI_, kind="linear", fill_value="extrapolate")
gNVII = interp1d(Temp, gNVII_, kind="linear", fill_value="extrapolate")
gNVIII = interp1d(Temp, gNVIII_, kind="linear", fill_value="extrapolate")
#-------------------------------------------


#----- cooling_rate_4d 
def cooling_rate_4d(ionX, T, nHI, nelec, nHII, Temp_4d, HIDensity_4d, elecDensity_4d, HIIDensity_4d):

  T = np.log10(T)
  nHI = np.log10(nHI)
  nelec = np.log10(nelec)
  nHII = np.log10(nHII)
  
  if ionX == 'CI':
    j = 0
  if ionX == 'FeI':
    j = 1

  if T <= 4:
    CI_rates = rates_4d[j, :]
    interp_4d = RegularGridInterpolator((Temp_4d, HIDensity_4d, elecDensity_4d, HIIDensity_4d), CI_rates)
    res = interp_4d(np.array([T, nHI, nelec, nHII]))[0]
    res = np.log10(10**res / 10**nelec) # because ne is multiplied later it will be cancelled by this. Intrinsically ne was included in chimes Table!!!
  else:
    CI_rates = rates_hiT_4d[j, :]
    interp_4d = interp1d(Temp_hiT_4d, CI_rates, kind='linear', fill_value="extrapolate")
    res = interp_4d(T)

  return res

#----- cooling_rate_2d 
def cooling_rate_2d(ionX, T, nelec, Temp_2d, elecDensity_2d): # include Temp_hiT here !!!!!!!!!!!!!!!!!!!!!!!

  T = np.log10(T)
  nelec = np.log10(nelec)
  
  if ionX == 'CII':
    j = 0
  if ionX == 'NII':
    j = 1
  if ionX == 'SiII':
    j = 2
  if ionX == 'FeII':
    j = 3

  if T <= 4.95:
    CII_rates = rates_2d[j, :]
    interp_2d = RegularGridInterpolator((Temp_2d, elecDensity_2d), CII_rates)
    res = interp_2d(np.array([T, nelec]))[0]
  else:
    CII_rates = rates_hiT_2d[j, :]
    interp_2d = interp1d(Temp_hiT_2d, CII_rates, kind='linear', fill_value="extrapolate")
    res = interp_2d(T)

  return res


#----- grain_recomb_rate
def grain_recomb_rate(ionx, T, nelec, G0, A_v, Temp, Psi): # Note: Temp and Psi are in log from CHIMES table!
  
  elmz = np.array(['HII', 'HeII', 'CII', 'OII', 'SiII', 'FeII', 'MgII', 'SII', 'CaII', 'CaIII'])
  ni = np.where(elmz == ionx)[0][0]
  
  Psix = G0 * np.exp(-2.77 * A_v) * T**0.5 / nelec
  Psix = np.log10(Psix)
  
  T = np.log10(T)
  
  interp_2d = RegularGridInterpolator((Temp, Psi), grain_recomb_rates[ni, :])
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
  
  

