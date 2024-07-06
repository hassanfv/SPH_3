
# In this version, I include grain_recombination (29 May 2024).

import h5py
from scipy.interpolate import interp1d
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import interp2d

gamma = 5./3.
kB = 1.3807e-16

with h5py.File('../chimes_main_data.hdf5', 'r') as file:
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



#-------------------------------------------
#------------ Carbon Section ---------------
#-------------------------------------------
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



#-------------------------------------------
#----------- Nitrogen Section --------------
#-------------------------------------------
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



#-----------------------------------------
#----------- Oxygen Section --------------
#-----------------------------------------
#--- REACTION RATES
R_OI_to_OII_via_HeII_ = rates[173, :]
R_OI_to_OII_via_e_ = rates[174, :]
R_OI_to_OII_via_HII_ = rates[181, :]
R_OII_to_OI_via_HI_ = rates[169, :]
R_OII_to_OI_via_e_ = rates[170, :]
R_OII_to_OIII_via_e_ = rates[171, :]
R_OIII_to_OII_via_HeI_ = rates[165, :]
R_OIII_to_OII_via_HI_ = rates[166, :]
R_OIII_to_OII_via_e_ = rates[167, :]
R_OIII_to_OIV_via_e_ = rates[168, :]
R_OIV_to_OV_via_e_ = rates[172, :]
R_OIV_to_OIII_via_e_ = rates[182, :]
R_OIV_to_OIII_via_HI_ = rates[183, :]
R_OIV_to_OIII_via_HeI_ = rates[184, :]
R_OV_to_OIV_via_HeI_ = rates[198, :]
R_OV_to_OIV_via_HI_ = rates[199, :]
R_OV_to_OIV_via_e_ = rates[200, :]
R_OV_to_OVI_via_e_ = rates[201, :]
R_OVI_to_OV_via_HI_ = rates[195, :]
R_OVI_to_OV_via_e_ = rates[196, :]
R_OVI_to_OVII_via_e_ = rates[197, :]
R_OVII_to_OVI_via_e_ = rates[193, :]
R_OVII_to_OVIII_via_e_ = rates[194, :]
R_OVIII_to_OVII_via_e_ = rates[191, :]
R_OVIII_to_OIX_via_e_ = rates[192, :]
R_OIX_to_OVIII_via_e_ = rates[190, :]
R_Om_to_OI_via_HII_ = rates[104, :]

R_OI_to_OII_via_HeII = interp1d(Temp, R_OI_to_OII_via_HeII_, kind="linear", fill_value="extrapolate")
R_OI_to_OII_via_e = interp1d(Temp, R_OI_to_OII_via_e_, kind="linear", fill_value="extrapolate")
R_OI_to_OII_via_HII = interp1d(Temp, R_OI_to_OII_via_HII_, kind="linear", fill_value="extrapolate")
R_OII_to_OI_via_HI = interp1d(Temp, R_OII_to_OI_via_HI_, kind="linear", fill_value="extrapolate")
R_OII_to_OI_via_e = interp1d(Temp, R_OII_to_OI_via_e_, kind="linear", fill_value="extrapolate")
R_OII_to_OIII_via_e = interp1d(Temp, R_OII_to_OIII_via_e_, kind="linear", fill_value="extrapolate")
R_OIII_to_OII_via_HeI = interp1d(Temp, R_OIII_to_OII_via_HeI_, kind="linear", fill_value="extrapolate")
R_OIII_to_OII_via_HI = interp1d(Temp, R_OIII_to_OII_via_HI_, kind="linear", fill_value="extrapolate")
R_OIII_to_OII_via_e = interp1d(Temp, R_OIII_to_OII_via_e_, kind="linear", fill_value="extrapolate")
R_OIII_to_OIV_via_e = interp1d(Temp, R_OIII_to_OIV_via_e_, kind="linear", fill_value="extrapolate")
R_OIV_to_OV_via_e = interp1d(Temp, R_OIV_to_OV_via_e_, kind="linear", fill_value="extrapolate")
R_OIV_to_OIII_via_e = interp1d(Temp, R_OIV_to_OIII_via_e_, kind="linear", fill_value="extrapolate")
R_OIV_to_OIII_via_HI = interp1d(Temp, R_OIV_to_OIII_via_HI_, kind="linear", fill_value="extrapolate")
R_OIV_to_OIII_via_HeI = interp1d(Temp, R_OIV_to_OIII_via_HeI_, kind="linear", fill_value="extrapolate")
R_OV_to_OIV_via_HeI = interp1d(Temp, R_OV_to_OIV_via_HeI_, kind="linear", fill_value="extrapolate")
R_OV_to_OIV_via_HI = interp1d(Temp, R_OV_to_OIV_via_HI_, kind="linear", fill_value="extrapolate")
R_OV_to_OIV_via_e = interp1d(Temp, R_OV_to_OIV_via_e_, kind="linear", fill_value="extrapolate")
R_OV_to_OVI_via_e = interp1d(Temp, R_OV_to_OVI_via_e_, kind="linear", fill_value="extrapolate")
R_OVI_to_OV_via_HI = interp1d(Temp, R_OVI_to_OV_via_HI_, kind="linear", fill_value="extrapolate")
R_OVI_to_OV_via_e = interp1d(Temp, R_OVI_to_OV_via_e_, kind="linear", fill_value="extrapolate")
R_OVI_to_OVII_via_e = interp1d(Temp, R_OVI_to_OVII_via_e_, kind="linear", fill_value="extrapolate")
R_OVII_to_OVI_via_e = interp1d(Temp, R_OVII_to_OVI_via_e_, kind="linear", fill_value="extrapolate")
R_OVII_to_OVIII_via_e = interp1d(Temp, R_OVII_to_OVIII_via_e_, kind="linear", fill_value="extrapolate")
R_OVIII_to_OVII_via_e = interp1d(Temp, R_OVIII_to_OVII_via_e_, kind="linear", fill_value="extrapolate")
R_OVIII_to_OIX_via_e = interp1d(Temp, R_OVIII_to_OIX_via_e_, kind="linear", fill_value="extrapolate")
R_OIX_to_OVIII_via_e = interp1d(Temp, R_OIX_to_OVIII_via_e_, kind="linear", fill_value="extrapolate")
R_Om_to_OI_via_HII = interp1d(Temp, R_Om_to_OI_via_HII_, kind="linear", fill_value="extrapolate")

#--- COOLING RATES
gOII_ = cooling_rates[17, :]
gOIII_ = cooling_rates[18, :]
gOIV_ = cooling_rates[19, :]
gOV_ = cooling_rates[20, :]
gOVI_ = cooling_rates[21, :]
gOVII_ = cooling_rates[22, :]
gOVIII_ = cooling_rates[23, :]
gOIX_ = cooling_rates[24, :]

gOII = interp1d(Temp, gOII_, kind="linear", fill_value="extrapolate")
gOIII = interp1d(Temp, gOIII_, kind="linear", fill_value="extrapolate")
gOIV = interp1d(Temp, gOIV_, kind="linear", fill_value="extrapolate")
gOV = interp1d(Temp, gOV_, kind="linear", fill_value="extrapolate")
gOVI = interp1d(Temp, gOVI_, kind="linear", fill_value="extrapolate")
gOVII = interp1d(Temp, gOVII_, kind="linear", fill_value="extrapolate")
gOVIII = interp1d(Temp, gOVIII_, kind="linear", fill_value="extrapolate")
gOIX = interp1d(Temp, gOIX_, kind="linear", fill_value="extrapolate")



#-----------------------------------------
#------------ Neon Section ---------------
#-----------------------------------------
#--- REACTION RATES
R_NeI_to_NeII_via_e_ = rates[189, :]
R_NeII_to_NeI_via_e_ = rates[187, :]
R_NeII_to_NeIII_via_e_ = rates[188, :]
R_NeIII_to_NeII_via_e_ = rates[185, :]
R_NeIII_to_NeIV_via_e_ = rates[186, :]
R_NeIII_to_NeII_via_HI_ = rates[224, :]
R_NeIII_to_NeII_via_HeI_ = rates[245, :]
R_NeIV_to_NeIII_via_e_ = rates[247, :]
R_NeIV_to_NeV_via_e_ = rates[269, :]
R_NeIV_to_NeIII_via_HeI_ = rates[302, :]
R_NeIV_to_NeIII_via_HI_ = rates[304, :]
R_NeV_to_NeIV_via_HeI_ = rates[298, :]
R_NeV_to_NeIV_via_HI_ = rates[299, :]
R_NeV_to_NeIV_via_e_ = rates[300, :]
R_NeV_to_NeVI_via_e_ = rates[301, :]
R_NeVI_to_NeV_via_HI_ = rates[295, :]
R_NeVI_to_NeV_via_e_ = rates[296, :]
R_NeVI_to_NeVII_via_e_ = rates[297, :]
R_NeVII_to_NeVI_via_e_ = rates[293, :]
R_NeVII_to_NeVIII_via_e_ = rates[294, :]
R_NeVIII_to_NeVII_via_e_ = rates[291, :]
R_NeVIII_to_NeIX_via_e_ = rates[292, :]
R_NeIX_to_NeVIII_via_e_ = rates[289, :]
R_NeIX_to_NeX_via_e_ = rates[290, :]
R_NeX_to_NeXI_via_e_ = rates[288, :]
R_NeX_to_NeIX_via_e_ = rates[305, :]
R_NeXI_to_NeX_via_e_ = rates[306, :]

R_NeI_to_NeII_via_e = interp1d(Temp, R_NeI_to_NeII_via_e_, kind="linear", fill_value="extrapolate")
R_NeII_to_NeI_via_e = interp1d(Temp, R_NeII_to_NeI_via_e_, kind="linear", fill_value="extrapolate")
R_NeII_to_NeIII_via_e = interp1d(Temp, R_NeII_to_NeIII_via_e_, kind="linear", fill_value="extrapolate")
R_NeIII_to_NeII_via_e = interp1d(Temp, R_NeIII_to_NeII_via_e_, kind="linear", fill_value="extrapolate")
R_NeIII_to_NeIV_via_e = interp1d(Temp, R_NeIII_to_NeIV_via_e_, kind="linear", fill_value="extrapolate")
R_NeIII_to_NeII_via_HI = interp1d(Temp, R_NeIII_to_NeII_via_HI_, kind="linear", fill_value="extrapolate")
R_NeIII_to_NeII_via_HeI = interp1d(Temp, R_NeIII_to_NeII_via_HeI_, kind="linear", fill_value="extrapolate")
R_NeIV_to_NeIII_via_e = interp1d(Temp, R_NeIV_to_NeIII_via_e_, kind="linear", fill_value="extrapolate")
R_NeIV_to_NeV_via_e = interp1d(Temp, R_NeIV_to_NeV_via_e_, kind="linear", fill_value="extrapolate")
R_NeIV_to_NeIII_via_HeI = interp1d(Temp, R_NeIV_to_NeIII_via_HeI_, kind="linear", fill_value="extrapolate")
R_NeIV_to_NeIII_via_HI = interp1d(Temp, R_NeIV_to_NeIII_via_HI_, kind="linear", fill_value="extrapolate")
R_NeV_to_NeIV_via_HeI = interp1d(Temp, R_NeV_to_NeIV_via_HeI_, kind="linear", fill_value="extrapolate")
R_NeV_to_NeIV_via_HI = interp1d(Temp, R_NeV_to_NeIV_via_HI_, kind="linear", fill_value="extrapolate")
R_NeV_to_NeIV_via_e = interp1d(Temp, R_NeV_to_NeIV_via_e_, kind="linear", fill_value="extrapolate")
R_NeV_to_NeVI_via_e = interp1d(Temp, R_NeV_to_NeVI_via_e_, kind="linear", fill_value="extrapolate")
R_NeVI_to_NeV_via_HI = interp1d(Temp, R_NeVI_to_NeV_via_HI_, kind="linear", fill_value="extrapolate")
R_NeVI_to_NeV_via_e = interp1d(Temp, R_NeVI_to_NeV_via_e_, kind="linear", fill_value="extrapolate")
R_NeVI_to_NeVII_via_e = interp1d(Temp, R_NeVI_to_NeVII_via_e_, kind="linear", fill_value="extrapolate")
R_NeVII_to_NeVI_via_e = interp1d(Temp, R_NeVII_to_NeVI_via_e_, kind="linear", fill_value="extrapolate")
R_NeVII_to_NeVIII_via_e = interp1d(Temp, R_NeVII_to_NeVIII_via_e_, kind="linear", fill_value="extrapolate")
R_NeVIII_to_NeVII_via_e = interp1d(Temp, R_NeVIII_to_NeVII_via_e_, kind="linear", fill_value="extrapolate")
R_NeVIII_to_NeIX_via_e = interp1d(Temp, R_NeVIII_to_NeIX_via_e_, kind="linear", fill_value="extrapolate")
R_NeIX_to_NeVIII_via_e = interp1d(Temp, R_NeIX_to_NeVIII_via_e_, kind="linear", fill_value="extrapolate")
R_NeIX_to_NeX_via_e = interp1d(Temp, R_NeIX_to_NeX_via_e_, kind="linear", fill_value="extrapolate")
R_NeX_to_NeXI_via_e = interp1d(Temp, R_NeX_to_NeXI_via_e_, kind="linear", fill_value="extrapolate")
R_NeX_to_NeIX_via_e = interp1d(Temp, R_NeX_to_NeIX_via_e_, kind="linear", fill_value="extrapolate")
R_NeXI_to_NeX_via_e = interp1d(Temp, R_NeXI_to_NeX_via_e_, kind="linear", fill_value="extrapolate")

#--- COOLING RATES
gNeI_ = cooling_rates[25, :]
gNeII_ = cooling_rates[26, :]
gNeIII_ = cooling_rates[27, :]
gNeIV_ = cooling_rates[28, :]
gNeV_ = cooling_rates[29, :]
gNeVI_ = cooling_rates[30, :]
gNeVII_ = cooling_rates[31, :]
gNeVIII_ = cooling_rates[32, :]
gNeIX_ = cooling_rates[33, :]
gNeX_ = cooling_rates[34, :]
gNeXI_ = cooling_rates[35, :]

gNeI = interp1d(Temp, gNeI_, kind="linear", fill_value="extrapolate")
gNeII = interp1d(Temp, gNeII_, kind="linear", fill_value="extrapolate")
gNeIII = interp1d(Temp, gNeIII_, kind="linear", fill_value="extrapolate")
gNeIV = interp1d(Temp, gNeIV_, kind="linear", fill_value="extrapolate")
gNeV = interp1d(Temp, gNeV_, kind="linear", fill_value="extrapolate")
gNeVI = interp1d(Temp, gNeVI_, kind="linear", fill_value="extrapolate")
gNeVII = interp1d(Temp, gNeVII_, kind="linear", fill_value="extrapolate")
gNeVIII = interp1d(Temp, gNeVIII_, kind="linear", fill_value="extrapolate")
gNeIX = interp1d(Temp, gNeIX_, kind="linear", fill_value="extrapolate")
gNeX = interp1d(Temp, gNeX_, kind="linear", fill_value="extrapolate")
gNeXI = interp1d(Temp, gNeXI_, kind="linear", fill_value="extrapolate")



#----------------------------------------------
#------------ Magnesium Section ---------------
#----------------------------------------------
#--- REACTION RATES
R_MgI_to_MgII_via_HII_ = rates[308, :]
R_MgI_to_MgII_via_e_ = rates[309, :]
R_MgII_to_MgIII_via_HII_ = rates[323, :]
R_MgII_to_MgI_via_e_ = rates[324, :]
R_MgII_to_MgIII_via_e_ = rates[326, :]
R_MgIII_to_MgII_via_HI_ = rates[320, :]
R_MgIII_to_MgII_via_e_ = rates[321, :]
R_MgIII_to_MgIV_via_e_ = rates[322, :]
R_MgIV_to_MgIII_via_HeI_ = rates[316, :]
R_MgIV_to_MgIII_via_e_ = rates[318, :]
R_MgIV_to_MgV_via_e_ = rates[319, :]
R_MgIV_to_MgIII_via_HI_ = rates[325, :]
R_MgV_to_MgIV_via_HeI_ = rates[312, :]
R_MgV_to_MgIV_via_HI_ = rates[313, :]
R_MgV_to_MgIV_via_e_ = rates[314, :]
R_MgV_to_MgVI_via_e_ = rates[315, :]
R_MgVI_to_MgV_via_HI_ = rates[310, :]
R_MgVI_to_MgVII_via_e_ = rates[311, :]
R_MgVI_to_MgV_via_e_ = rates[317, :]
R_MgVII_to_MgVIII_via_e_ = rates[287, :]
R_MgVII_to_MgVI_via_e_ = rates[303, :]
R_MgVIII_to_MgVII_via_e_ = rates[265, :]
R_MgVIII_to_MgIX_via_e_ = rates[285, :]
R_MgIX_to_MgVIII_via_e_ = rates[262, :]
R_MgIX_to_MgX_via_e_ = rates[263, :]
R_MgX_to_MgIX_via_e_ = rates[260, :]
R_MgX_to_MgXI_via_e_ = rates[261, :]
R_MgXI_to_MgX_via_e_ = rates[258, :]
R_MgXI_to_MgXII_via_e_ = rates[259, :]
R_MgXII_to_MgXIII_via_e_ = rates[257, :]
R_MgXII_to_MgXI_via_e_ = rates[264, :]
R_MgXIII_to_MgXII_via_e_ = rates[256, :]

R_MgI_to_MgII_via_HII = interp1d(Temp, R_MgI_to_MgII_via_HII_, kind="linear", fill_value="extrapolate")
R_MgI_to_MgII_via_e = interp1d(Temp, R_MgI_to_MgII_via_e_, kind="linear", fill_value="extrapolate")
R_MgII_to_MgIII_via_HII = interp1d(Temp, R_MgII_to_MgIII_via_HII_, kind="linear", fill_value="extrapolate")
R_MgII_to_MgI_via_e = interp1d(Temp, R_MgII_to_MgI_via_e_, kind="linear", fill_value="extrapolate")
R_MgII_to_MgIII_via_e = interp1d(Temp, R_MgII_to_MgIII_via_e_, kind="linear", fill_value="extrapolate")
R_MgIII_to_MgII_via_HI = interp1d(Temp, R_MgIII_to_MgII_via_HI_, kind="linear", fill_value="extrapolate")
R_MgIII_to_MgII_via_e = interp1d(Temp, R_MgIII_to_MgII_via_e_, kind="linear", fill_value="extrapolate")
R_MgIII_to_MgIV_via_e = interp1d(Temp, R_MgIII_to_MgIV_via_e_, kind="linear", fill_value="extrapolate")
R_MgIV_to_MgIII_via_HeI = interp1d(Temp, R_MgIV_to_MgIII_via_HeI_, kind="linear", fill_value="extrapolate")
R_MgIV_to_MgIII_via_e = interp1d(Temp, R_MgIV_to_MgIII_via_e_, kind="linear", fill_value="extrapolate")
R_MgIV_to_MgV_via_e = interp1d(Temp, R_MgIV_to_MgV_via_e_, kind="linear", fill_value="extrapolate")
R_MgIV_to_MgIII_via_HI = interp1d(Temp, R_MgIV_to_MgIII_via_HI_, kind="linear", fill_value="extrapolate")
R_MgV_to_MgIV_via_HeI = interp1d(Temp, R_MgV_to_MgIV_via_HeI_, kind="linear", fill_value="extrapolate")
R_MgV_to_MgIV_via_HI = interp1d(Temp, R_MgV_to_MgIV_via_HI_, kind="linear", fill_value="extrapolate")
R_MgV_to_MgIV_via_e = interp1d(Temp, R_MgV_to_MgIV_via_e_, kind="linear", fill_value="extrapolate")
R_MgV_to_MgVI_via_e = interp1d(Temp, R_MgV_to_MgVI_via_e_, kind="linear", fill_value="extrapolate")
R_MgVI_to_MgV_via_HI = interp1d(Temp, R_MgVI_to_MgV_via_HI_, kind="linear", fill_value="extrapolate")
R_MgVI_to_MgVII_via_e = interp1d(Temp, R_MgVI_to_MgVII_via_e_, kind="linear", fill_value="extrapolate")
R_MgVI_to_MgV_via_e = interp1d(Temp, R_MgVI_to_MgV_via_e_, kind="linear", fill_value="extrapolate")
R_MgVII_to_MgVIII_via_e = interp1d(Temp, R_MgVII_to_MgVIII_via_e_, kind="linear", fill_value="extrapolate")
R_MgVII_to_MgVI_via_e = interp1d(Temp, R_MgVII_to_MgVI_via_e_, kind="linear", fill_value="extrapolate")
R_MgVIII_to_MgVII_via_e = interp1d(Temp, R_MgVIII_to_MgVII_via_e_, kind="linear", fill_value="extrapolate")
R_MgVIII_to_MgIX_via_e = interp1d(Temp, R_MgVIII_to_MgIX_via_e_, kind="linear", fill_value="extrapolate")
R_MgIX_to_MgVIII_via_e = interp1d(Temp, R_MgIX_to_MgVIII_via_e_, kind="linear", fill_value="extrapolate")
R_MgIX_to_MgX_via_e = interp1d(Temp, R_MgIX_to_MgX_via_e_, kind="linear", fill_value="extrapolate")
R_MgX_to_MgIX_via_e = interp1d(Temp, R_MgX_to_MgIX_via_e_, kind="linear", fill_value="extrapolate")
R_MgX_to_MgXI_via_e = interp1d(Temp, R_MgX_to_MgXI_via_e_, kind="linear", fill_value="extrapolate")
R_MgXI_to_MgX_via_e = interp1d(Temp, R_MgXI_to_MgX_via_e_, kind="linear", fill_value="extrapolate")
R_MgXI_to_MgXII_via_e = interp1d(Temp, R_MgXI_to_MgXII_via_e_, kind="linear", fill_value="extrapolate")
R_MgXII_to_MgXIII_via_e = interp1d(Temp, R_MgXII_to_MgXIII_via_e_, kind="linear", fill_value="extrapolate")
R_MgXII_to_MgXI_via_e = interp1d(Temp, R_MgXII_to_MgXI_via_e_, kind="linear", fill_value="extrapolate")
R_MgXIII_to_MgXII_via_e = interp1d(Temp, R_MgXIII_to_MgXII_via_e_, kind="linear", fill_value="extrapolate")

#--- COOLING RATES
gMgI_ = cooling_rates[36, :]
gMgII_ = cooling_rates[37, :]
gMgIII_ = cooling_rates[38, :]
gMgIV_ = cooling_rates[39, :]
gMgV_ = cooling_rates[40, :]
gMgVI_ = cooling_rates[41, :]
gMgVII_ = cooling_rates[42, :]
gMgVIII_ = cooling_rates[43, :]
gMgIX_ = cooling_rates[44, :]
gMgX_ = cooling_rates[45, :]
gMgXI_ = cooling_rates[46, :]
gMgXII_ = cooling_rates[47, :]
gMgXIII_ = cooling_rates[48, :]

gMgI = interp1d(Temp, gMgI_, kind="linear", fill_value="extrapolate")
gMgII = interp1d(Temp, gMgII_, kind="linear", fill_value="extrapolate")
gMgIII = interp1d(Temp, gMgIII_, kind="linear", fill_value="extrapolate")
gMgIV = interp1d(Temp, gMgIV_, kind="linear", fill_value="extrapolate")
gMgV = interp1d(Temp, gMgV_, kind="linear", fill_value="extrapolate")
gMgVI = interp1d(Temp, gMgVI_, kind="linear", fill_value="extrapolate")
gMgVII = interp1d(Temp, gMgVII_, kind="linear", fill_value="extrapolate")
gMgVIII = interp1d(Temp, gMgVIII_, kind="linear", fill_value="extrapolate")
gMgIX = interp1d(Temp, gMgIX_, kind="linear", fill_value="extrapolate")
gMgX = interp1d(Temp, gMgX_, kind="linear", fill_value="extrapolate")
gMgXI = interp1d(Temp, gMgXI_, kind="linear", fill_value="extrapolate")
gMgXII = interp1d(Temp, gMgXII_, kind="linear", fill_value="extrapolate")
gMgXIII = interp1d(Temp, gMgXIII_, kind="linear", fill_value="extrapolate")



#--------------------------------------------
#------------ Silicon Section ---------------
#--------------------------------------------
#--- REACTION RATES
R_SiI_to_SiII_via_HeII_ = rates[252, :]
R_SiI_to_SiII_via_HII_ = rates[253, :]
R_SiI_to_SiII_via_e_ = rates[254, :]
R_SiII_to_SiIII_via_HII_ = rates[249, :]
R_SiII_to_SiI_via_e_ = rates[250, :]
R_SiII_to_SiIII_via_e_ = rates[251, :]
R_SiII_to_SiIII_via_HeII_ = rates[286, :]
R_SiIII_to_SiIV_via_e_ = rates[255, :]
R_SiIII_to_SiII_via_e_ = rates[266, :]
R_SiIII_to_SiII_via_HI_ = rates[267, :]
R_SiIII_to_SiIV_via_HeII_ = rates[275, :]
R_SiIV_to_SiIII_via_HeI_ = rates[280, :]
R_SiIV_to_SiIII_via_HI_ = rates[281, :]
R_SiIV_to_SiIII_via_e_ = rates[282, :]
R_SiIV_to_SiV_via_e_ = rates[283, :]
R_SiV_to_SiIV_via_HI_ = rates[277, :]
R_SiV_to_SiIV_via_e_ = rates[278, :]
R_SiV_to_SiVI_via_e_ = rates[279, :]
R_SiV_to_SiIV_via_HeI_ = rates[284, :]
R_SiVI_to_SiV_via_HI_ = rates[273, :]
R_SiVI_to_SiV_via_e_ = rates[274, :]
R_SiVI_to_SiVII_via_e_ = rates[276, :]
R_SiVII_to_SiVI_via_e_ = rates[271, :]
R_SiVII_to_SiVIII_via_e_ = rates[272, :]
R_SiVIII_to_SiVII_via_e_ = rates[268, :]
R_SiVIII_to_SiIX_via_e_ = rates[270, :]
R_SiIX_to_SiX_via_e_ = rates[248, :]
R_SiIX_to_SiVIII_via_e_ = rates[307, :]
R_SiX_to_SiXI_via_e_ = rates[164, :]
R_SiX_to_SiIX_via_e_ = rates[179, :]
R_SiXI_to_SiX_via_e_ = rates[80, :]
R_SiXI_to_SiXII_via_e_ = rates[162, :]
R_SiXII_to_SiXI_via_e_ = rates[56, :]
R_SiXII_to_SiXIII_via_e_ = rates[57, :]
R_SiXIII_to_SiXII_via_e_ = rates[54, :]
R_SiXIII_to_SiXIV_via_e_ = rates[55, :]
R_SiXIV_to_SiXIII_via_e_ = rates[52, :]
R_SiXIV_to_SiXV_via_e_ = rates[53, :]
R_SiXV_to_SiXIV_via_e_ = rates[51, :]

R_SiI_to_SiII_via_HeII = interp1d(Temp, R_SiI_to_SiII_via_HeII_, kind="linear", fill_value="extrapolate")
R_SiI_to_SiII_via_HII = interp1d(Temp, R_SiI_to_SiII_via_HII_, kind="linear", fill_value="extrapolate")
R_SiI_to_SiII_via_e = interp1d(Temp, R_SiI_to_SiII_via_e_, kind="linear", fill_value="extrapolate")
R_SiII_to_SiIII_via_HII = interp1d(Temp, R_SiII_to_SiIII_via_HII_, kind="linear", fill_value="extrapolate")
R_SiII_to_SiI_via_e = interp1d(Temp, R_SiII_to_SiI_via_e_, kind="linear", fill_value="extrapolate")
R_SiII_to_SiIII_via_e = interp1d(Temp, R_SiII_to_SiIII_via_e_, kind="linear", fill_value="extrapolate")
R_SiII_to_SiIII_via_HeII = interp1d(Temp, R_SiII_to_SiIII_via_HeII_, kind="linear", fill_value="extrapolate")
R_SiIII_to_SiIV_via_e = interp1d(Temp, R_SiIII_to_SiIV_via_e_, kind="linear", fill_value="extrapolate")
R_SiIII_to_SiII_via_e = interp1d(Temp, R_SiIII_to_SiII_via_e_, kind="linear", fill_value="extrapolate")
R_SiIII_to_SiII_via_HI = interp1d(Temp, R_SiIII_to_SiII_via_HI_, kind="linear", fill_value="extrapolate")
R_SiIII_to_SiIV_via_HeII = interp1d(Temp, R_SiIII_to_SiIV_via_HeII_, kind="linear", fill_value="extrapolate")
R_SiIV_to_SiIII_via_HeI = interp1d(Temp, R_SiIV_to_SiIII_via_HeI_, kind="linear", fill_value="extrapolate")
R_SiIV_to_SiIII_via_HI = interp1d(Temp, R_SiIV_to_SiIII_via_HI_, kind="linear", fill_value="extrapolate")
R_SiIV_to_SiIII_via_e = interp1d(Temp, R_SiIV_to_SiIII_via_e_, kind="linear", fill_value="extrapolate")
R_SiIV_to_SiV_via_e = interp1d(Temp, R_SiIV_to_SiV_via_e_, kind="linear", fill_value="extrapolate")
R_SiV_to_SiIV_via_HI = interp1d(Temp, R_SiV_to_SiIV_via_HI_, kind="linear", fill_value="extrapolate")
R_SiV_to_SiIV_via_e = interp1d(Temp, R_SiV_to_SiIV_via_e_, kind="linear", fill_value="extrapolate")
R_SiV_to_SiVI_via_e = interp1d(Temp, R_SiV_to_SiVI_via_e_, kind="linear", fill_value="extrapolate")
R_SiV_to_SiIV_via_HeI = interp1d(Temp, R_SiV_to_SiIV_via_HeI_, kind="linear", fill_value="extrapolate")
R_SiVI_to_SiV_via_HI = interp1d(Temp, R_SiVI_to_SiV_via_HI_, kind="linear", fill_value="extrapolate")
R_SiVI_to_SiV_via_e = interp1d(Temp, R_SiVI_to_SiV_via_e_, kind="linear", fill_value="extrapolate")
R_SiVI_to_SiVII_via_e = interp1d(Temp, R_SiVI_to_SiVII_via_e_, kind="linear", fill_value="extrapolate")
R_SiVII_to_SiVI_via_e = interp1d(Temp, R_SiVII_to_SiVI_via_e_, kind="linear", fill_value="extrapolate")
R_SiVII_to_SiVIII_via_e = interp1d(Temp, R_SiVII_to_SiVIII_via_e_, kind="linear", fill_value="extrapolate")
R_SiVIII_to_SiVII_via_e = interp1d(Temp, R_SiVIII_to_SiVII_via_e_, kind="linear", fill_value="extrapolate")
R_SiVIII_to_SiIX_via_e = interp1d(Temp, R_SiVIII_to_SiIX_via_e_, kind="linear", fill_value="extrapolate")
R_SiIX_to_SiX_via_e = interp1d(Temp, R_SiIX_to_SiX_via_e_, kind="linear", fill_value="extrapolate")
R_SiIX_to_SiVIII_via_e = interp1d(Temp, R_SiIX_to_SiVIII_via_e_, kind="linear", fill_value="extrapolate")
R_SiX_to_SiXI_via_e = interp1d(Temp, R_SiX_to_SiXI_via_e_, kind="linear", fill_value="extrapolate")
R_SiX_to_SiIX_via_e = interp1d(Temp, R_SiX_to_SiIX_via_e_, kind="linear", fill_value="extrapolate")
R_SiXI_to_SiX_via_e = interp1d(Temp, R_SiXI_to_SiX_via_e_, kind="linear", fill_value="extrapolate")
R_SiXI_to_SiXII_via_e = interp1d(Temp, R_SiXI_to_SiXII_via_e_, kind="linear", fill_value="extrapolate")
R_SiXII_to_SiXI_via_e = interp1d(Temp, R_SiXII_to_SiXI_via_e_, kind="linear", fill_value="extrapolate")
R_SiXII_to_SiXIII_via_e = interp1d(Temp, R_SiXII_to_SiXIII_via_e_, kind="linear", fill_value="extrapolate")
R_SiXIII_to_SiXII_via_e = interp1d(Temp, R_SiXIII_to_SiXII_via_e_, kind="linear", fill_value="extrapolate")
R_SiXIII_to_SiXIV_via_e = interp1d(Temp, R_SiXIII_to_SiXIV_via_e_, kind="linear", fill_value="extrapolate")
R_SiXIV_to_SiXIII_via_e = interp1d(Temp, R_SiXIV_to_SiXIII_via_e_, kind="linear", fill_value="extrapolate")
R_SiXIV_to_SiXV_via_e = interp1d(Temp, R_SiXIV_to_SiXV_via_e_, kind="linear", fill_value="extrapolate")
R_SiXV_to_SiXIV_via_e = interp1d(Temp, R_SiXV_to_SiXIV_via_e_, kind="linear", fill_value="extrapolate")

#--- COOLING RATES
gSiI_ = cooling_rates[49, :]
gSiIII_ = cooling_rates[50, :]
gSiIV_ = cooling_rates[51, :]
gSiV_ = cooling_rates[52, :]
gSiVI_ = cooling_rates[53, :]
gSiVII_ = cooling_rates[54, :]
gSiVIII_ = cooling_rates[55, :]
gSiIX_ = cooling_rates[56, :]
gSiX_ = cooling_rates[57, :]
gSiXI_ = cooling_rates[58, :]
gSiXII_ = cooling_rates[59, :]
gSiXIII_ = cooling_rates[60, :]
gSiXIV_ = cooling_rates[61, :]
gSiXV_ = cooling_rates[62, :]

gSiI = interp1d(Temp, gSiI_, kind="linear", fill_value="extrapolate")
gSiIII = interp1d(Temp, gSiIII_, kind="linear", fill_value="extrapolate")
gSiIV = interp1d(Temp, gSiIV_, kind="linear", fill_value="extrapolate")
gSiV = interp1d(Temp, gSiV_, kind="linear", fill_value="extrapolate")
gSiVI = interp1d(Temp, gSiVI_, kind="linear", fill_value="extrapolate")
gSiVII = interp1d(Temp, gSiVII_, kind="linear", fill_value="extrapolate")
gSiVIII = interp1d(Temp, gSiVIII_, kind="linear", fill_value="extrapolate")
gSiIX = interp1d(Temp, gSiIX_, kind="linear", fill_value="extrapolate")
gSiX = interp1d(Temp, gSiX_, kind="linear", fill_value="extrapolate")
gSiXI = interp1d(Temp, gSiXI_, kind="linear", fill_value="extrapolate")
gSiXII = interp1d(Temp, gSiXII_, kind="linear", fill_value="extrapolate")
gSiXIII = interp1d(Temp, gSiXIII_, kind="linear", fill_value="extrapolate")
gSiXIV = interp1d(Temp, gSiXIV_, kind="linear", fill_value="extrapolate")
gSiXV = interp1d(Temp, gSiXV_, kind="linear", fill_value="extrapolate")



#-----------------------------------------
#------------ Sulpher Section ---------------
#-----------------------------------------
#--- REACTION RATES
R_SI_to_SII_via_HII_ = rates[50, :]
R_SI_to_SII_via_e_ = rates[58, :]
R_SII_to_SI_via_HI_ = rates[45, :]
R_SII_to_SIII_via_HeII_ = rates[46, :]
R_SII_to_SI_via_e_ = rates[47, :]
R_SII_to_SIII_via_e_ = rates[48, :]
R_SIII_to_SIV_via_HeII_ = rates[42, :]
R_SIII_to_SII_via_e_ = rates[43, :]
R_SIII_to_SIV_via_e_ = rates[44, :]
R_SIII_to_SII_via_HI_ = rates[49, :]
R_SIV_to_SV_via_e_ = rates[59, :]
R_SIV_to_SIII_via_e_ = rates[60, :]
R_SIV_to_SIII_via_HI_ = rates[61, :]
R_SIV_to_SIII_via_HeI_ = rates[78, :]
R_SV_to_SIV_via_HeI_ = rates[74, :]
R_SV_to_SIV_via_HI_ = rates[75, :]
R_SV_to_SIV_via_e_ = rates[76, :]
R_SV_to_SVI_via_e_ = rates[77, :]
R_SVI_to_SV_via_HI_ = rates[71, :]
R_SVI_to_SV_via_e_ = rates[72, :]
R_SVI_to_SVII_via_e_ = rates[73, :]
R_SVII_to_SVI_via_e_ = rates[69, :]
R_SVII_to_SVIII_via_e_ = rates[70, :]
R_SVIII_to_SVII_via_e_ = rates[67, :]
R_SVIII_to_SIX_via_e_ = rates[68, :]
R_SIX_to_SVIII_via_e_ = rates[65, :]
R_SIX_to_SX_via_e_ = rates[66, :]
R_SX_to_SXI_via_e_ = rates[64, :]
R_SX_to_SIX_via_e_ = rates[163, :]
R_SXI_to_SX_via_e_ = rates[62, :]
R_SXI_to_SXII_via_e_ = rates[63, :]
R_SXII_to_SXIII_via_e_ = rates[41, :]
R_SXII_to_SXI_via_e_ = rates[79, :]
R_SXIII_to_SXII_via_e_ = rates[38, :]
R_SXIII_to_SXIV_via_e_ = rates[40, :]
R_SXIV_to_SXIII_via_e_ = rates[15, :]
R_SXIV_to_SXV_via_e_ = rates[16, :]
R_SXV_to_SXIV_via_e_ = rates[13, :]
R_SXV_to_SXVI_via_e_ = rates[14, :]
R_SXVI_to_SXV_via_e_ = rates[11, :]
R_SXVI_to_SXVII_via_e_ = rates[12, :]
R_SXVII_to_SXVI_via_e_ = rates[10, :]


R_SI_to_SII_via_HII = interp1d(Temp, R_SI_to_SII_via_HII_, kind="linear", fill_value="extrapolate")
R_SI_to_SII_via_e = interp1d(Temp, R_SI_to_SII_via_e_, kind="linear", fill_value="extrapolate")
R_SII_to_SI_via_HI = interp1d(Temp, R_SII_to_SI_via_HI_, kind="linear", fill_value="extrapolate")
R_SII_to_SIII_via_HeII = interp1d(Temp, R_SII_to_SIII_via_HeII_, kind="linear", fill_value="extrapolate")
R_SII_to_SI_via_e = interp1d(Temp, R_SII_to_SI_via_e_, kind="linear", fill_value="extrapolate")
R_SII_to_SIII_via_e = interp1d(Temp, R_SII_to_SIII_via_e_, kind="linear", fill_value="extrapolate")
R_SIII_to_SIV_via_HeII = interp1d(Temp, R_SIII_to_SIV_via_HeII_, kind="linear", fill_value="extrapolate")
R_SIII_to_SII_via_e = interp1d(Temp, R_SIII_to_SII_via_e_, kind="linear", fill_value="extrapolate")
R_SIII_to_SIV_via_e = interp1d(Temp, R_SIII_to_SIV_via_e_, kind="linear", fill_value="extrapolate")
R_SIII_to_SII_via_HI = interp1d(Temp, R_SIII_to_SII_via_HI_, kind="linear", fill_value="extrapolate")
R_SIV_to_SV_via_e = interp1d(Temp, R_SIV_to_SV_via_e_, kind="linear", fill_value="extrapolate")
R_SIV_to_SIII_via_e = interp1d(Temp, R_SIV_to_SIII_via_e_, kind="linear", fill_value="extrapolate")
R_SIV_to_SIII_via_HI = interp1d(Temp, R_SIV_to_SIII_via_HI_, kind="linear", fill_value="extrapolate")
R_SIV_to_SIII_via_HeI = interp1d(Temp, R_SIV_to_SIII_via_HeI_, kind="linear", fill_value="extrapolate")
R_SV_to_SIV_via_HeI = interp1d(Temp, R_SV_to_SIV_via_HeI_, kind="linear", fill_value="extrapolate")
R_SV_to_SIV_via_HI = interp1d(Temp, R_SV_to_SIV_via_HI_, kind="linear", fill_value="extrapolate")
R_SV_to_SIV_via_e = interp1d(Temp, R_SV_to_SIV_via_e_, kind="linear", fill_value="extrapolate")
R_SV_to_SVI_via_e = interp1d(Temp, R_SV_to_SVI_via_e_, kind="linear", fill_value="extrapolate")
R_SVI_to_SV_via_HI = interp1d(Temp, R_SVI_to_SV_via_HI_, kind="linear", fill_value="extrapolate")
R_SVI_to_SV_via_e = interp1d(Temp, R_SVI_to_SV_via_e_, kind="linear", fill_value="extrapolate")
R_SVI_to_SVII_via_e = interp1d(Temp, R_SVI_to_SVII_via_e_, kind="linear", fill_value="extrapolate")
R_SVII_to_SVI_via_e = interp1d(Temp, R_SVII_to_SVI_via_e_, kind="linear", fill_value="extrapolate")
R_SVII_to_SVIII_via_e = interp1d(Temp, R_SVII_to_SVIII_via_e_, kind="linear", fill_value="extrapolate")
R_SVIII_to_SVII_via_e = interp1d(Temp, R_SVIII_to_SVII_via_e_, kind="linear", fill_value="extrapolate")
R_SVIII_to_SIX_via_e = interp1d(Temp, R_SVIII_to_SIX_via_e_, kind="linear", fill_value="extrapolate")
R_SIX_to_SVIII_via_e = interp1d(Temp, R_SIX_to_SVIII_via_e_, kind="linear", fill_value="extrapolate")
R_SIX_to_SX_via_e = interp1d(Temp, R_SIX_to_SX_via_e_, kind="linear", fill_value="extrapolate")
R_SX_to_SXI_via_e = interp1d(Temp, R_SX_to_SXI_via_e_, kind="linear", fill_value="extrapolate")
R_SX_to_SIX_via_e = interp1d(Temp, R_SX_to_SIX_via_e_, kind="linear", fill_value="extrapolate")
R_SXI_to_SX_via_e = interp1d(Temp, R_SXI_to_SX_via_e_, kind="linear", fill_value="extrapolate")
R_SXI_to_SXII_via_e = interp1d(Temp, R_SXI_to_SXII_via_e_, kind="linear", fill_value="extrapolate")
R_SXII_to_SXIII_via_e = interp1d(Temp, R_SXII_to_SXIII_via_e_, kind="linear", fill_value="extrapolate")
R_SXII_to_SXI_via_e = interp1d(Temp, R_SXII_to_SXI_via_e_, kind="linear", fill_value="extrapolate")
R_SXIII_to_SXII_via_e = interp1d(Temp, R_SXIII_to_SXII_via_e_, kind="linear", fill_value="extrapolate")
R_SXIII_to_SXIV_via_e = interp1d(Temp, R_SXIII_to_SXIV_via_e_, kind="linear", fill_value="extrapolate")
R_SXIV_to_SXIII_via_e = interp1d(Temp, R_SXIV_to_SXIII_via_e_, kind="linear", fill_value="extrapolate")
R_SXIV_to_SXV_via_e = interp1d(Temp, R_SXIV_to_SXV_via_e_, kind="linear", fill_value="extrapolate")
R_SXV_to_SXIV_via_e = interp1d(Temp, R_SXV_to_SXIV_via_e_, kind="linear", fill_value="extrapolate")
R_SXV_to_SXVI_via_e = interp1d(Temp, R_SXV_to_SXVI_via_e_, kind="linear", fill_value="extrapolate")
R_SXVI_to_SXV_via_e = interp1d(Temp, R_SXVI_to_SXV_via_e_, kind="linear", fill_value="extrapolate")
R_SXVI_to_SXVII_via_e = interp1d(Temp, R_SXVI_to_SXVII_via_e_, kind="linear", fill_value="extrapolate")
R_SXVII_to_SXVI_via_e = interp1d(Temp, R_SXVII_to_SXVI_via_e_, kind="linear", fill_value="extrapolate")

#--- COOLING RATES
gSI_ = cooling_rates[63, :]
gSII_ = cooling_rates[64, :]
gSIII_ = cooling_rates[65, :]
gSIV_ = cooling_rates[66, :]
gSV_ = cooling_rates[67, :]
gSVI_ = cooling_rates[68, :]
gSVII_ = cooling_rates[69, :]
gSVIII_ = cooling_rates[70, :]
gSIX_ = cooling_rates[71, :]
gSX_ = cooling_rates[72, :]
gSXI_ = cooling_rates[73, :]
gSXII_ = cooling_rates[74, :]
gSXIII_ = cooling_rates[75, :]
gSXIV_ = cooling_rates[76, :]
gSXV_ = cooling_rates[77, :]
gSXVI_ = cooling_rates[78, :]
gSXVII_ = cooling_rates[79, :]

gSI = interp1d(Temp, gSI_, kind="linear", fill_value="extrapolate")
gSII = interp1d(Temp, gSII_, kind="linear", fill_value="extrapolate")
gSIII = interp1d(Temp, gSIII_, kind="linear", fill_value="extrapolate")
gSIV = interp1d(Temp, gSIV_, kind="linear", fill_value="extrapolate")
gSV = interp1d(Temp, gSV_, kind="linear", fill_value="extrapolate")
gSVI = interp1d(Temp, gSVI_, kind="linear", fill_value="extrapolate")
gSVII = interp1d(Temp, gSVII_, kind="linear", fill_value="extrapolate")
gSVIII = interp1d(Temp, gSVIII_, kind="linear", fill_value="extrapolate")
gSIX = interp1d(Temp, gSIX_, kind="linear", fill_value="extrapolate")
gSX = interp1d(Temp, gSX_, kind="linear", fill_value="extrapolate")
gSXI = interp1d(Temp, gSXI_, kind="linear", fill_value="extrapolate")
gSXII = interp1d(Temp, gSXII_, kind="linear", fill_value="extrapolate")
gSXIII = interp1d(Temp, gSXIII_, kind="linear", fill_value="extrapolate")
gSXIV = interp1d(Temp, gSXIV_, kind="linear", fill_value="extrapolate")
gSXV = interp1d(Temp, gSXV_, kind="linear", fill_value="extrapolate")
gSXVI = interp1d(Temp, gSXVI_, kind="linear", fill_value="extrapolate")
gSXVII = interp1d(Temp, gSXVII_, kind="linear", fill_value="extrapolate")


#--------------------------------------------
#------------ Calcium Section ---------------
#--------------------------------------------
#--- REACTION RATES
R_CaI_to_CaII_via_e_ = rates[17, :]
R_CaII_to_CaI_via_e_ = rates[7, :]
R_CaII_to_CaIII_via_e_ = rates[9, :]
R_CaIII_to_CaII_via_e_ = rates[5, :]
R_CaIII_to_CaIV_via_e_ = rates[6, :]
R_CaIV_to_CaIII_via_HI_ = rates[2, :]
R_CaIV_to_CaIII_via_e_ = rates[3, :]
R_CaIV_to_CaV_via_e_ = rates[4, :]
R_CaV_to_CaVI_via_e_ = rates[1, :]
R_CaV_to_CaIV_via_e_ = rates[8, :]
R_CaV_to_CaIV_via_HI_ = rates[18, :]
R_CaVI_to_CaVII_via_e_ = rates[19, :]
R_CaVI_to_CaV_via_e_ = rates[20, :]
R_CaVI_to_CaV_via_HI_ = rates[37, :]
R_CaVII_to_CaVI_via_e_ = rates[35, :]
R_CaVII_to_CaVIII_via_e_ = rates[36, :]
R_CaVIII_to_CaVII_via_e_ = rates[33, :]
R_CaVIII_to_CaIX_via_e_ = rates[34, :]
R_CaIX_to_CaVIII_via_e_ = rates[31, :]
R_CaIX_to_CaX_via_e_ = rates[32, :]
R_CaX_to_CaIX_via_e_ = rates[29, :]
R_CaX_to_CaXI_via_e_ = rates[30, :]
R_CaXI_to_CaX_via_e_ = rates[27, :]
R_CaXI_to_CaXII_via_e_ = rates[28, :]
R_CaXII_to_CaXI_via_e_ = rates[25, :]
R_CaXII_to_CaXIII_via_e_ = rates[26, :]
R_CaXIII_to_CaXII_via_e_ = rates[23, :]
R_CaXIII_to_CaXIV_via_e_ = rates[24, :]
R_CaXIV_to_CaXIII_via_e_ = rates[21, :]
R_CaXIV_to_CaXV_via_e_ = rates[22, :]
R_CaXV_to_CaXVI_via_e_ = rates[39, :]
R_CaXV_to_CaXIV_via_e_ = rates[81, :]
R_CaXVI_to_CaXV_via_e_ = rates[82, :]
R_CaXVI_to_CaXVII_via_e_ = rates[121, :]
R_CaXVII_to_CaXVI_via_e_ = rates[138, :]
R_CaXVII_to_CaXVIII_via_e_ = rates[139, :]
R_CaXVIII_to_CaXVII_via_e_ = rates[136, :]
R_CaXVIII_to_CaXIX_via_e_ = rates[137, :]
R_CaXIX_to_CaXVIII_via_e_ = rates[134, :]
R_CaXIX_to_CaXX_via_e_ = rates[135, :]
R_CaXX_to_CaXXI_via_e_ = rates[133, :]
R_CaXX_to_CaXIX_via_e_ = rates[140, :]
R_CaXXI_to_CaXX_via_e_ = rates[132, :]

R_CaI_to_CaII_via_e = interp1d(Temp, R_CaI_to_CaII_via_e_, kind="linear", fill_value="extrapolate")
R_CaII_to_CaI_via_e = interp1d(Temp, R_CaII_to_CaI_via_e_, kind="linear", fill_value="extrapolate")
R_CaII_to_CaIII_via_e = interp1d(Temp, R_CaII_to_CaIII_via_e_, kind="linear", fill_value="extrapolate")
R_CaIII_to_CaII_via_e = interp1d(Temp, R_CaIII_to_CaII_via_e_, kind="linear", fill_value="extrapolate")
R_CaIII_to_CaIV_via_e = interp1d(Temp, R_CaIII_to_CaIV_via_e_, kind="linear", fill_value="extrapolate")
R_CaIV_to_CaIII_via_HI = interp1d(Temp, R_CaIV_to_CaIII_via_HI_, kind="linear", fill_value="extrapolate")
R_CaIV_to_CaIII_via_e = interp1d(Temp, R_CaIV_to_CaIII_via_e_, kind="linear", fill_value="extrapolate")
R_CaIV_to_CaV_via_e = interp1d(Temp, R_CaIV_to_CaV_via_e_, kind="linear", fill_value="extrapolate")
R_CaV_to_CaVI_via_e = interp1d(Temp, R_CaV_to_CaVI_via_e_, kind="linear", fill_value="extrapolate")
R_CaV_to_CaIV_via_e = interp1d(Temp, R_CaV_to_CaIV_via_e_, kind="linear", fill_value="extrapolate")
R_CaV_to_CaIV_via_HI = interp1d(Temp, R_CaV_to_CaIV_via_HI_, kind="linear", fill_value="extrapolate")
R_CaVI_to_CaVII_via_e = interp1d(Temp, R_CaVI_to_CaVII_via_e_, kind="linear", fill_value="extrapolate")
R_CaVI_to_CaV_via_e = interp1d(Temp, R_CaVI_to_CaV_via_e_, kind="linear", fill_value="extrapolate")
R_CaVI_to_CaV_via_HI = interp1d(Temp, R_CaVI_to_CaV_via_HI_, kind="linear", fill_value="extrapolate")
R_CaVII_to_CaVI_via_e = interp1d(Temp, R_CaVII_to_CaVI_via_e_, kind="linear", fill_value="extrapolate")
R_CaVII_to_CaVIII_via_e = interp1d(Temp, R_CaVII_to_CaVIII_via_e_, kind="linear", fill_value="extrapolate")
R_CaVIII_to_CaVII_via_e = interp1d(Temp, R_CaVIII_to_CaVII_via_e_, kind="linear", fill_value="extrapolate")
R_CaVIII_to_CaIX_via_e = interp1d(Temp, R_CaVIII_to_CaIX_via_e_, kind="linear", fill_value="extrapolate")
R_CaIX_to_CaVIII_via_e = interp1d(Temp, R_CaIX_to_CaVIII_via_e_, kind="linear", fill_value="extrapolate")
R_CaIX_to_CaX_via_e = interp1d(Temp, R_CaIX_to_CaX_via_e_, kind="linear", fill_value="extrapolate")
R_CaX_to_CaIX_via_e = interp1d(Temp, R_CaX_to_CaIX_via_e_, kind="linear", fill_value="extrapolate")
R_CaX_to_CaXI_via_e = interp1d(Temp, R_CaX_to_CaXI_via_e_, kind="linear", fill_value="extrapolate")
R_CaXI_to_CaX_via_e = interp1d(Temp, R_CaXI_to_CaX_via_e_, kind="linear", fill_value="extrapolate")
R_CaXI_to_CaXII_via_e = interp1d(Temp, R_CaXI_to_CaXII_via_e_, kind="linear", fill_value="extrapolate")
R_CaXII_to_CaXI_via_e = interp1d(Temp, R_CaXII_to_CaXI_via_e_, kind="linear", fill_value="extrapolate")
R_CaXII_to_CaXIII_via_e = interp1d(Temp, R_CaXII_to_CaXIII_via_e_, kind="linear", fill_value="extrapolate")
R_CaXIII_to_CaXII_via_e = interp1d(Temp, R_CaXIII_to_CaXII_via_e_, kind="linear", fill_value="extrapolate")
R_CaXIII_to_CaXIV_via_e = interp1d(Temp, R_CaXIII_to_CaXIV_via_e_, kind="linear", fill_value="extrapolate")
R_CaXIV_to_CaXIII_via_e = interp1d(Temp, R_CaXIV_to_CaXIII_via_e_, kind="linear", fill_value="extrapolate")
R_CaXIV_to_CaXV_via_e = interp1d(Temp, R_CaXIV_to_CaXV_via_e_, kind="linear", fill_value="extrapolate")
R_CaXV_to_CaXVI_via_e = interp1d(Temp, R_CaXV_to_CaXVI_via_e_, kind="linear", fill_value="extrapolate")
R_CaXV_to_CaXIV_via_e = interp1d(Temp, R_CaXV_to_CaXIV_via_e_, kind="linear", fill_value="extrapolate")
R_CaXVI_to_CaXV_via_e = interp1d(Temp, R_CaXVI_to_CaXV_via_e_, kind="linear", fill_value="extrapolate")
R_CaXVI_to_CaXVII_via_e = interp1d(Temp, R_CaXVI_to_CaXVII_via_e_, kind="linear", fill_value="extrapolate")
R_CaXVII_to_CaXVI_via_e = interp1d(Temp, R_CaXVII_to_CaXVI_via_e_, kind="linear", fill_value="extrapolate")
R_CaXVII_to_CaXVIII_via_e = interp1d(Temp, R_CaXVII_to_CaXVIII_via_e_, kind="linear", fill_value="extrapolate")
R_CaXVIII_to_CaXVII_via_e = interp1d(Temp, R_CaXVIII_to_CaXVII_via_e_, kind="linear", fill_value="extrapolate")
R_CaXVIII_to_CaXIX_via_e = interp1d(Temp, R_CaXVIII_to_CaXIX_via_e_, kind="linear", fill_value="extrapolate")
R_CaXIX_to_CaXVIII_via_e = interp1d(Temp, R_CaXIX_to_CaXVIII_via_e_, kind="linear", fill_value="extrapolate")
R_CaXIX_to_CaXX_via_e = interp1d(Temp, R_CaXIX_to_CaXX_via_e_, kind="linear", fill_value="extrapolate")
R_CaXX_to_CaXXI_via_e = interp1d(Temp, R_CaXX_to_CaXXI_via_e_, kind="linear", fill_value="extrapolate")
R_CaXX_to_CaXIX_via_e = interp1d(Temp, R_CaXX_to_CaXIX_via_e_, kind="linear", fill_value="extrapolate")
R_CaXXI_to_CaXX_via_e = interp1d(Temp, R_CaXXI_to_CaXX_via_e_, kind="linear", fill_value="extrapolate")

#--- COOLING RATES
gCaI_ = cooling_rates[80, :]
gCaII_ = cooling_rates[81, :]
gCaIII_ = cooling_rates[82, :]
gCaIV_ = cooling_rates[83, :]
gCaV_ = cooling_rates[84, :]
gCaVI_ = cooling_rates[85, :]
gCaVII_ = cooling_rates[86, :]
gCaVIII_ = cooling_rates[87, :]
gCaIX_ = cooling_rates[88, :]
gCaX_ = cooling_rates[89, :]
gCaXI_ = cooling_rates[90, :]
gCaXII_ = cooling_rates[91, :]
gCaXIII_ = cooling_rates[92, :]
gCaXIV_ = cooling_rates[93, :]
gCaXV_ = cooling_rates[94, :]
gCaXVI_ = cooling_rates[95, :]
gCaXVII_ = cooling_rates[96, :]
gCaXVIII_ = cooling_rates[97, :]
gCaXIX_ = cooling_rates[98, :]
gCaXX_ = cooling_rates[99, :]
gCaXXI_ = cooling_rates[100, :]

gCaI = interp1d(Temp, gCaI_, kind="linear", fill_value="extrapolate")
gCaII = interp1d(Temp, gCaII_, kind="linear", fill_value="extrapolate")
gCaIII = interp1d(Temp, gCaIII_, kind="linear", fill_value="extrapolate")
gCaIV = interp1d(Temp, gCaIV_, kind="linear", fill_value="extrapolate")
gCaV = interp1d(Temp, gCaV_, kind="linear", fill_value="extrapolate")
gCaVI = interp1d(Temp, gCaVI_, kind="linear", fill_value="extrapolate")
gCaVII = interp1d(Temp, gCaVII_, kind="linear", fill_value="extrapolate")
gCaVIII = interp1d(Temp, gCaVIII_, kind="linear", fill_value="extrapolate")
gCaIX = interp1d(Temp, gCaIX_, kind="linear", fill_value="extrapolate")
gCaX = interp1d(Temp, gCaX_, kind="linear", fill_value="extrapolate")
gCaXI = interp1d(Temp, gCaXI_, kind="linear", fill_value="extrapolate")
gCaXII = interp1d(Temp, gCaXII_, kind="linear", fill_value="extrapolate")
gCaXIII = interp1d(Temp, gCaXIII_, kind="linear", fill_value="extrapolate")
gCaXIV = interp1d(Temp, gCaXIV_, kind="linear", fill_value="extrapolate")
gCaXV = interp1d(Temp, gCaXV_, kind="linear", fill_value="extrapolate")
gCaXVI = interp1d(Temp, gCaXVI_, kind="linear", fill_value="extrapolate")
gCaXVII = interp1d(Temp, gCaXVII_, kind="linear", fill_value="extrapolate")
gCaXVIII = interp1d(Temp, gCaXVIII_, kind="linear", fill_value="extrapolate")
gCaXIX = interp1d(Temp, gCaXIX_, kind="linear", fill_value="extrapolate")
gCaXX = interp1d(Temp, gCaXX_, kind="linear", fill_value="extrapolate")
gCaXXI = interp1d(Temp, gCaXXI_, kind="linear", fill_value="extrapolate")



#--------------------------------------------
#-------------- Iron Section ----------------
#--------------------------------------------
#--- REACTION RATES
R_FeI_to_FeII_via_HII_ = rates[129, :]
R_FeI_to_FeII_via_e_ = rates[130, :]
R_FeII_to_FeIII_via_HII_ = rates[126, :]
R_FeII_to_FeI_via_e_ = rates[127, :]
R_FeII_to_FeIII_via_e_ = rates[128, :]
R_FeIII_to_FeII_via_e_ = rates[124, :]
R_FeIII_to_FeIV_via_e_ = rates[125, :]
R_FeIII_to_FeII_via_HI_ = rates[131, :]
R_FeIV_to_FeV_via_e_ = rates[141, :]
R_FeIV_to_FeIII_via_e_ = rates[142, :]
R_FeIV_to_FeIII_via_HI_ = rates[143, :]
R_FeV_to_FeIV_via_HI_ = rates[158, :]
R_FeV_to_FeIV_via_e_ = rates[159, :]
R_FeV_to_FeVI_via_e_ = rates[160, :]
R_FeVI_to_FeV_via_HI_ = rates[155, :]
R_FeVI_to_FeV_via_e_ = rates[156, :]
R_FeVI_to_FeVII_via_e_ = rates[157, :]
R_FeVII_to_FeVI_via_e_ = rates[153, :]
R_FeVII_to_FeVIII_via_e_ = rates[154, :]
R_FeVIII_to_FeVII_via_e_ = rates[151, :]
R_FeVIII_to_FeIX_via_e_ = rates[152, :]
R_FeIX_to_FeVIII_via_e_ = rates[149, :]
R_FeIX_to_FeX_via_e_ = rates[150, :]
R_FeX_to_FeIX_via_e_ = rates[147, :]
R_FeX_to_FeXI_via_e_ = rates[148, :]
R_FeXI_to_FeX_via_e_ = rates[145, :]
R_FeXI_to_FeXII_via_e_ = rates[146, :]
R_FeXII_to_FeXI_via_e_ = rates[123, :]
R_FeXII_to_FeXIII_via_e_ = rates[144, :]
R_FeXIII_to_FeXII_via_e_ = rates[122, :]
R_FeXIII_to_FeXIV_via_e_ = rates[161, :]
R_FeXIV_to_FeXIII_via_e_ = rates[98, :]
R_FeXIV_to_FeXV_via_e_ = rates[120, :]
R_FeXV_to_FeXIV_via_e_ = rates[96, :]
R_FeXV_to_FeXVI_via_e_ = rates[97, :]
R_FeXVI_to_FeXV_via_e_ = rates[94, :]
R_FeXVI_to_FeXVII_via_e_ = rates[95, :]
R_FeXVII_to_FeXVI_via_e_ = rates[92, :]
R_FeXVII_to_FeXVIII_via_e_ = rates[93, :]
R_FeXVIII_to_FeXVII_via_e_ = rates[91, :]
R_FeXVIII_to_FeXIX_via_e_ = rates[99, :]
R_FeXIX_to_FeXVIII_via_e_ = rates[88, :]
R_FeXIX_to_FeXX_via_e_ = rates[89, :]
R_FeXX_to_FeXIX_via_e_ = rates[86, :]
R_FeXX_to_FeXXI_via_e_ = rates[87, :]
R_FeXXI_to_FeXX_via_e_ = rates[84, :]
R_FeXXI_to_FeXXII_via_e_ = rates[85, :]
R_FeXXII_to_FeXXIII_via_e_ = rates[83, :]
R_FeXXII_to_FeXXI_via_e_ = rates[90, :]
R_FeXXIII_to_FeXXIV_via_e_ = rates[100, :]
R_FeXXIII_to_FeXXII_via_e_ = rates[101, :]
R_FeXXIV_to_FeXXV_via_e_ = rates[102, :]
R_FeXXIV_to_FeXXIII_via_e_ = rates[119, :]
R_FeXXV_to_FeXXIV_via_e_ = rates[117, :]
R_FeXXV_to_FeXXVI_via_e_ = rates[118, :]
R_FeXXVI_to_FeXXV_via_e_ = rates[115, :]
R_FeXXVI_to_FeXXVII_via_e_ = rates[116, :]
R_FeXXVII_to_FeXXVI_via_e_ = rates[114, :]

R_FeI_to_FeII_via_HII = interp1d(Temp, R_FeI_to_FeII_via_HII_, kind="linear", fill_value="extrapolate")
R_FeI_to_FeII_via_e = interp1d(Temp, R_FeI_to_FeII_via_e_, kind="linear", fill_value="extrapolate")
R_FeII_to_FeIII_via_HII = interp1d(Temp, R_FeII_to_FeIII_via_HII_, kind="linear", fill_value="extrapolate")
R_FeII_to_FeI_via_e = interp1d(Temp, R_FeII_to_FeI_via_e_, kind="linear", fill_value="extrapolate")
R_FeII_to_FeIII_via_e = interp1d(Temp, R_FeII_to_FeIII_via_e_, kind="linear", fill_value="extrapolate")
R_FeIII_to_FeII_via_e = interp1d(Temp, R_FeIII_to_FeII_via_e_, kind="linear", fill_value="extrapolate")
R_FeIII_to_FeIV_via_e = interp1d(Temp, R_FeIII_to_FeIV_via_e_, kind="linear", fill_value="extrapolate")
R_FeIII_to_FeII_via_HI = interp1d(Temp, R_FeIII_to_FeII_via_HI_, kind="linear", fill_value="extrapolate")
R_FeIV_to_FeV_via_e = interp1d(Temp, R_FeIV_to_FeV_via_e_, kind="linear", fill_value="extrapolate")
R_FeIV_to_FeIII_via_e = interp1d(Temp, R_FeIV_to_FeIII_via_e_, kind="linear", fill_value="extrapolate")
R_FeIV_to_FeIII_via_HI = interp1d(Temp, R_FeIV_to_FeIII_via_HI_, kind="linear", fill_value="extrapolate")
R_FeV_to_FeIV_via_HI = interp1d(Temp, R_FeV_to_FeIV_via_HI_, kind="linear", fill_value="extrapolate")
R_FeV_to_FeIV_via_e = interp1d(Temp, R_FeV_to_FeIV_via_e_, kind="linear", fill_value="extrapolate")
R_FeV_to_FeVI_via_e = interp1d(Temp, R_FeV_to_FeVI_via_e_, kind="linear", fill_value="extrapolate")
R_FeVI_to_FeV_via_HI = interp1d(Temp, R_FeVI_to_FeV_via_HI_, kind="linear", fill_value="extrapolate")
R_FeVI_to_FeV_via_e = interp1d(Temp, R_FeVI_to_FeV_via_e_, kind="linear", fill_value="extrapolate")
R_FeVI_to_FeVII_via_e = interp1d(Temp, R_FeVI_to_FeVII_via_e_, kind="linear", fill_value="extrapolate")
R_FeVII_to_FeVI_via_e = interp1d(Temp, R_FeVII_to_FeVI_via_e_, kind="linear", fill_value="extrapolate")
R_FeVII_to_FeVIII_via_e = interp1d(Temp, R_FeVII_to_FeVIII_via_e_, kind="linear", fill_value="extrapolate")
R_FeVIII_to_FeVII_via_e = interp1d(Temp, R_FeVIII_to_FeVII_via_e_, kind="linear", fill_value="extrapolate")
R_FeVIII_to_FeIX_via_e = interp1d(Temp, R_FeVIII_to_FeIX_via_e_, kind="linear", fill_value="extrapolate")
R_FeIX_to_FeVIII_via_e = interp1d(Temp, R_FeIX_to_FeVIII_via_e_, kind="linear", fill_value="extrapolate")
R_FeIX_to_FeX_via_e = interp1d(Temp, R_FeIX_to_FeX_via_e_, kind="linear", fill_value="extrapolate")
R_FeX_to_FeIX_via_e = interp1d(Temp, R_FeX_to_FeIX_via_e_, kind="linear", fill_value="extrapolate")
R_FeX_to_FeXI_via_e = interp1d(Temp, R_FeX_to_FeXI_via_e_, kind="linear", fill_value="extrapolate")
R_FeXI_to_FeX_via_e = interp1d(Temp, R_FeXI_to_FeX_via_e_, kind="linear", fill_value="extrapolate")
R_FeXI_to_FeXII_via_e = interp1d(Temp, R_FeXI_to_FeXII_via_e_, kind="linear", fill_value="extrapolate")
R_FeXII_to_FeXI_via_e = interp1d(Temp, R_FeXII_to_FeXI_via_e_, kind="linear", fill_value="extrapolate")
R_FeXII_to_FeXIII_via_e = interp1d(Temp, R_FeXII_to_FeXIII_via_e_, kind="linear", fill_value="extrapolate")
R_FeXIII_to_FeXII_via_e = interp1d(Temp, R_FeXIII_to_FeXII_via_e_, kind="linear", fill_value="extrapolate")
R_FeXIII_to_FeXIV_via_e = interp1d(Temp, R_FeXIII_to_FeXIV_via_e_, kind="linear", fill_value="extrapolate")
R_FeXIV_to_FeXIII_via_e = interp1d(Temp, R_FeXIV_to_FeXIII_via_e_, kind="linear", fill_value="extrapolate")
R_FeXIV_to_FeXV_via_e = interp1d(Temp, R_FeXIV_to_FeXV_via_e_, kind="linear", fill_value="extrapolate")
R_FeXV_to_FeXIV_via_e = interp1d(Temp, R_FeXV_to_FeXIV_via_e_, kind="linear", fill_value="extrapolate")
R_FeXV_to_FeXVI_via_e = interp1d(Temp, R_FeXV_to_FeXVI_via_e_, kind="linear", fill_value="extrapolate")
R_FeXVI_to_FeXV_via_e = interp1d(Temp, R_FeXVI_to_FeXV_via_e_, kind="linear", fill_value="extrapolate")
R_FeXVI_to_FeXVII_via_e = interp1d(Temp, R_FeXVI_to_FeXVII_via_e_, kind="linear", fill_value="extrapolate")
R_FeXVII_to_FeXVI_via_e = interp1d(Temp, R_FeXVII_to_FeXVI_via_e_, kind="linear", fill_value="extrapolate")
R_FeXVII_to_FeXVIII_via_e = interp1d(Temp, R_FeXVII_to_FeXVIII_via_e_, kind="linear", fill_value="extrapolate")
R_FeXVIII_to_FeXVII_via_e = interp1d(Temp, R_FeXVIII_to_FeXVII_via_e_, kind="linear", fill_value="extrapolate")
R_FeXVIII_to_FeXIX_via_e = interp1d(Temp, R_FeXVIII_to_FeXIX_via_e_, kind="linear", fill_value="extrapolate")
R_FeXIX_to_FeXVIII_via_e = interp1d(Temp, R_FeXIX_to_FeXVIII_via_e_, kind="linear", fill_value="extrapolate")
R_FeXIX_to_FeXX_via_e = interp1d(Temp, R_FeXIX_to_FeXX_via_e_, kind="linear", fill_value="extrapolate")
R_FeXX_to_FeXIX_via_e = interp1d(Temp, R_FeXX_to_FeXIX_via_e_, kind="linear", fill_value="extrapolate")
R_FeXX_to_FeXXI_via_e = interp1d(Temp, R_FeXX_to_FeXXI_via_e_, kind="linear", fill_value="extrapolate")
R_FeXXI_to_FeXX_via_e = interp1d(Temp, R_FeXXI_to_FeXX_via_e_, kind="linear", fill_value="extrapolate")
R_FeXXI_to_FeXXII_via_e = interp1d(Temp, R_FeXXI_to_FeXXII_via_e_, kind="linear", fill_value="extrapolate")
R_FeXXII_to_FeXXIII_via_e = interp1d(Temp, R_FeXXII_to_FeXXIII_via_e_, kind="linear", fill_value="extrapolate")
R_FeXXII_to_FeXXI_via_e = interp1d(Temp, R_FeXXII_to_FeXXI_via_e_, kind="linear", fill_value="extrapolate")
R_FeXXIII_to_FeXXIV_via_e = interp1d(Temp, R_FeXXIII_to_FeXXIV_via_e_, kind="linear", fill_value="extrapolate")
R_FeXXIII_to_FeXXII_via_e = interp1d(Temp, R_FeXXIII_to_FeXXII_via_e_, kind="linear", fill_value="extrapolate")
R_FeXXIV_to_FeXXV_via_e = interp1d(Temp, R_FeXXIV_to_FeXXV_via_e_, kind="linear", fill_value="extrapolate")
R_FeXXIV_to_FeXXIII_via_e = interp1d(Temp, R_FeXXIV_to_FeXXIII_via_e_, kind="linear", fill_value="extrapolate")
R_FeXXV_to_FeXXIV_via_e = interp1d(Temp, R_FeXXV_to_FeXXIV_via_e_, kind="linear", fill_value="extrapolate")
R_FeXXV_to_FeXXVI_via_e = interp1d(Temp, R_FeXXV_to_FeXXVI_via_e_, kind="linear", fill_value="extrapolate")
R_FeXXVI_to_FeXXV_via_e = interp1d(Temp, R_FeXXVI_to_FeXXV_via_e_, kind="linear", fill_value="extrapolate")
R_FeXXVI_to_FeXXVII_via_e = interp1d(Temp, R_FeXXVI_to_FeXXVII_via_e_, kind="linear", fill_value="extrapolate")
R_FeXXVII_to_FeXXVI_via_e = interp1d(Temp, R_FeXXVII_to_FeXXVI_via_e_, kind="linear", fill_value="extrapolate")

#--- COOLING RATES
gFeI_ = cooling_rates[101, :]
gFeIII_ = cooling_rates[102, :]
gFeIV_ = cooling_rates[103, :]
gFeV_ = cooling_rates[104, :]
gFeVI_ = cooling_rates[105, :]
gFeVII_ = cooling_rates[106, :]
gFeVIII_ = cooling_rates[107, :]
gFeIX_ = cooling_rates[108, :]
gFeX_ = cooling_rates[109, :]
gFeXI_ = cooling_rates[110, :]
gFeXII_ = cooling_rates[111, :]
gFeXIII_ = cooling_rates[112, :]
gFeXIV_ = cooling_rates[113, :]
gFeXV_ = cooling_rates[114, :]
gFeXVI_ = cooling_rates[115, :]
gFeXVII_ = cooling_rates[116, :]
gFeXVIII_ = cooling_rates[117, :]
gFeXIX_ = cooling_rates[118, :]
gFeXX_ = cooling_rates[119, :]
gFeXXI_ = cooling_rates[120, :]
gFeXXII_ = cooling_rates[121, :]
gFeXXIII_ = cooling_rates[122, :]
gFeXXIV_ = cooling_rates[123, :]
gFeXXV_ = cooling_rates[124, :]
gFeXXVI_ = cooling_rates[125, :]
gFeXXVII_ = cooling_rates[126, :]

gFeI = interp1d(Temp, gFeI_, kind="linear", fill_value="extrapolate")
gFeIII = interp1d(Temp, gFeIII_, kind="linear", fill_value="extrapolate")
gFeIV = interp1d(Temp, gFeIV_, kind="linear", fill_value="extrapolate")
gFeV = interp1d(Temp, gFeV_, kind="linear", fill_value="extrapolate")
gFeVI = interp1d(Temp, gFeVI_, kind="linear", fill_value="extrapolate")
gFeVII = interp1d(Temp, gFeVII_, kind="linear", fill_value="extrapolate")
gFeVIII = interp1d(Temp, gFeVIII_, kind="linear", fill_value="extrapolate")
gFeIX = interp1d(Temp, gFeIX_, kind="linear", fill_value="extrapolate")
gFeX = interp1d(Temp, gFeX_, kind="linear", fill_value="extrapolate")
gFeXI = interp1d(Temp, gFeXI_, kind="linear", fill_value="extrapolate")
gFeXII = interp1d(Temp, gFeXII_, kind="linear", fill_value="extrapolate")
gFeXIII = interp1d(Temp, gFeXIII_, kind="linear", fill_value="extrapolate")
gFeXIV = interp1d(Temp, gFeXIV_, kind="linear", fill_value="extrapolate")
gFeXV = interp1d(Temp, gFeXV_, kind="linear", fill_value="extrapolate")
gFeXVI = interp1d(Temp, gFeXVI_, kind="linear", fill_value="extrapolate")
gFeXVII = interp1d(Temp, gFeXVII_, kind="linear", fill_value="extrapolate")
gFeXVIII = interp1d(Temp, gFeXVIII_, kind="linear", fill_value="extrapolate")
gFeXIX = interp1d(Temp, gFeXIX_, kind="linear", fill_value="extrapolate")
gFeXX = interp1d(Temp, gFeXX_, kind="linear", fill_value="extrapolate")
gFeXXI = interp1d(Temp, gFeXXI_, kind="linear", fill_value="extrapolate")
gFeXXII = interp1d(Temp, gFeXXII_, kind="linear", fill_value="extrapolate")
gFeXXIII = interp1d(Temp, gFeXXIII_, kind="linear", fill_value="extrapolate")
gFeXXIV = interp1d(Temp, gFeXXIV_, kind="linear", fill_value="extrapolate")
gFeXXV = interp1d(Temp, gFeXXV_, kind="linear", fill_value="extrapolate")
gFeXXVI = interp1d(Temp, gFeXXVI_, kind="linear", fill_value="extrapolate")
gFeXXVII = interp1d(Temp, gFeXXVII_, kind="linear", fill_value="extrapolate")




#----- cooling_rate_4d 
def cooling_rate_4d(ionX, T, nHI, nelec, nHII, Temp_4d, HIDensity_4d, elecDensity_4d, HIIDensity_4d):

  T = np.log10(T)
  nHI = np.log10(nHI)
  nelec = np.log10(nelec)
  nHII = np.log10(nHII)
  
  if ionX == 'CI':
    j = 0
  if ionX == 'OI':
    j = 1

  if T <= 3.96:
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
  
  

