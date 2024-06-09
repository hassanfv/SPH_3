

# In this version, i.e. coolAuto_v5.py, the "func" function is automatically created by "create_func_script.py" code.
# In this version, i.e. coolAuto_v4.py, I included Auto-generated ODEs for H and He!
# In this version, i.e. coolAuto_v3.py, I defined the more general form of cooling_rate_4d and cooling_rate_2d functions!
# In this version, i.e. coolAuto_v2.py, I add Hm contribution.
# In this version, i.e. coolAuto_v1.py, the ODEs in func are created by the create_ode_script3_grain_const.py script!

import numpy as np
import matplotlib.pyplot as plt
#from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp
import pickle
from datax_v2 import *
from scipy.interpolate import RegularGridInterpolator


# Cooling due to free-free emission ::: will be multiplied by (nHp+nHep+nHepp)*ne later in the code.
def gfree(T):

  gff = 1.1 + 0.34 * np.exp(-(5.5 - np.log10(T))**2/3.0)
  g4_val = 1.42e-27 * gff * T**0.5
  
  return g4_val



#----- Lambda
def Lambda(T, nHI, nHII, nHm, nHeI, nHeII, nHeIII, nCI, nCII, nCIII, nCIV, nCV, 
           nCVI, nCVII, nNI, nNII, nNIII, nNIV, nNV, nNVI, nNVII, nNVIII, nCm):

  Tx = np.log10(T)

  ne = (
         1 * nHII - nHm + (nHeII + 2.0 * nHeIII) + 1 * nCII + 2 * nCIII + 3 * nCIV
       + 4 * nCV + 5 * nCVI + 6 * nCVII + 1 * nNII + 2 * nNIII
       + 3 * nNIV + 4 * nNV + 5 * nNVI + 6 * nNVII + 7 * nNVIII - 1 * nCm
       )

  cFree = (
            1 * nHII + nHeII + 4.0 * nHeIII + 1 * nCII + 4 * nCIII + 9 * nCIV
          + 16 * nCV + 25 * nCVI + 36 * nCVII + 1 * nNII + 4 * nNIII
          + 9 * nNIV + 16 * nNV + 25 * nNVI + 36 * nNVII + 49 * nNVIII )

  #----- # Glover & Jappsen - 2007 -----
  z = 0.0 # current time redshift!
  TCMB_0 = 2.7255
  TCMB = TCMB_0 * (1.0 + z)
  LCompton = 1.017e-37 * TCMB**4 * (T - TCMB) * ne
  #--------------------------------------

  grain_recomb_cooling_rate = grain_cool_rate(T, ne, nHI, nHII, G0, A_v, dust_ratio, Temp, Psi)

  Lamb = (
          10**g1(Tx) * ne * nHI  # H0
        + 10**g2(Tx) * ne * nHII # Hp
        + 10**grain_recomb_cooling_rate * nHII * ne # grain_recombination cooling!
        + 10**g3(Tx) * nHeI * ne # He0
        + 10**g4(Tx) * nHeII * ne # Hep
        + 10**grain_recomb_cooling_rate * nHeII * ne # grain_recombination cooling!
        + 10**g5(Tx) * nHeIII * ne# Hep
        + 10**cooling_rate_4d("CI", T, nHI, ne, nHII, Temp_4d, HIDensity_4d, elecDensity_4d, HIIDensity_4d) * nCI * ne # cooling via CI
        + 10**cooling_rate_2d("CII", T, ne, Temp_2d, elecDensity_2d) * nCII * ne # cooling via CII
        + 10**cooling_rate_2d("NII", T, ne, Temp_2d, elecDensity_2d) * nNII * ne # cooling via NII
        + 10**grain_recomb_cooling_rate * nCII * ne # grain_recombination cooling!
        + 10**gCIII(Tx) * nCIII * ne
        + 10**gCIV(Tx) * nCIV * ne
        + 10**gCV(Tx) * nCV * ne
        + 10**gCVI(Tx) * nCVI * ne
        + 10**gCVII(Tx) * nCVII * ne
        + 10**gNI(Tx) * nNI * ne
        + 10**gNIII(Tx) * nNIII * ne
        + 10**gNIV(Tx) * nNIV * ne
        + 10**gNV(Tx) * nNV * ne
        + 10**gNVI(Tx) * nNVI * ne
        + 10**gNVII(Tx) * nNVII * ne
        + 10**gNVIII(Tx) * nNVIII * ne
        + gfree(T) * ne * cFree # free-free emission
        + LCompton)

  return Lamb



#----- func
def func(t, y):

  nHI, nHII, nHm, nHeI, nHeII, nHeIII, nCI, nCII, nCIII, nCIV, \
  nCV, nCVI, nCVII, nCm, nNI, nNII, nNIII, nNIV, nNV, nNVI, \
  nNVII, nNVIII, T = y

  Tx = np.log10(T)

  ne = (
       + 1 * nHII + -1 * nHm + 1 * nHeII + 2 * nHeIII + 1 * nCII + 2 * nCIII
       + 3 * nCIV + 4 * nCV + 5 * nCVI + 6 * nCVII + -1 * nCm + 1 * nNII
       + 2 * nNIII + 3 * nNIV + 4 * nNV + 5 * nNVI + 6 * nNVII + 7 * nNVIII
      
       )

  ntot = (
           ne + nHI + nHII + nHm + nHeI + nHeII + nHeIII
         + nCI + nCII + nCIII + nCIV + nCV + nCVI
         + nCVII + nCm + nNI + nNII + nNIII + nNIV
         + nNV + nNVI + nNVII + nNVIII
         )

  grain_rec_HII_to_HI = grain_recomb_rate("HII", T, ne, G0, A_v, Temp, Psi)
  grain_rec_HeII_to_HeI = grain_recomb_rate("HeII", T, ne, G0, A_v, Temp, Psi)
  grain_rec_CII_to_CI = grain_recomb_rate("CII", T, ne, G0, A_v, Temp, Psi)

  const_CI_e_to_Cm_ = 2.2500E-15 # constant/rates

  dnHI_dt = (
             - 10**R_HI_to_HII_via_e(Tx) * nHI * ne
             + 10**grain_rec_HII_to_HI * nHII * ne # grain_recombination
             - 10**R_HI_to_Hm_via_e(Tx) * nHI * ne
             - 10**R_Hm_to_HI_via_HI(Tx) * nHm * nHI
             + 10**R_Hm_to_HI_via_HI(Tx) * nHm * nHI
             + 10**R_Hm_to_HI_via_e(Tx) * nHm * ne
             + 10**R_Hm_to_HI_via_HII(Tx) * nHm * nHII
             + 10**R_HeI_to_HeII_via_HII(Tx) * nHeI * nHII
             + 10**R_HeII_to_HeI_via_Hm(Tx) * nHeII * nHm
             - 10**R_HeII_to_HeI_via_HI(Tx) * nHeII * nHI
             - 10**R_HeIII_to_HeII_via_HI(Tx) * nHeIII * nHI
             + 10**R_CI_to_CII_via_HII(Tx) * nCI * nHII
             - 10**R_CII_to_CI_via_HI(Tx) * nCII * nHI
             - 10**R_CIII_to_CII_via_HI(Tx) * nCIII * nHI
             - 10**R_CIV_to_CIII_via_HI(Tx) * nCIV * nHI
             - 10**R_CV_to_CIV_via_HI(Tx) * nCV * nHI
             - 10**R_CVI_to_CV_via_HI(Tx) * nCVI * nHI
             + 10**R_Cm_to_CI_via_HII(Tx) * nCm * nHII
             + 10**R_NI_to_NII_via_HII(Tx) * nNI * nHII
             - 10**R_NII_to_NI_via_HI(Tx) * nNII * nHI
             - 10**R_NIII_to_NII_via_HI(Tx) * nNIII * nHI
             - 10**R_NIV_to_NIII_via_HI(Tx) * nNIV * nHI
             - 10**R_NV_to_NIV_via_HI(Tx) * nNV * nHI
             - 10**R_NVI_to_NV_via_HI(Tx) * nNVI * nHI
             + 10**R_HII_to_HI_via_e_caseA(Tx) * nHII * ne # H CaseA
            )

  dnHII_dt = (
              + 10**R_HI_to_HII_via_e(Tx) * nHI * ne
              - 10**grain_rec_HII_to_HI * nHII * ne # grain_recombination
              - 10**R_Hm_to_HI_via_HII(Tx) * nHm * nHII
              - 10**R_HeI_to_HeII_via_HII(Tx) * nHeI * nHII
              + 10**R_HeII_to_HeI_via_HI(Tx) * nHeII * nHI
              + 10**R_HeIII_to_HeII_via_HI(Tx) * nHeIII * nHI
              - 10**R_CI_to_CII_via_HII(Tx) * nCI * nHII
              + 10**R_CII_to_CI_via_HI(Tx) * nCII * nHI
              + 10**R_CIII_to_CII_via_HI(Tx) * nCIII * nHI
              + 10**R_CIV_to_CIII_via_HI(Tx) * nCIV * nHI
              + 10**R_CV_to_CIV_via_HI(Tx) * nCV * nHI
              + 10**R_CVI_to_CV_via_HI(Tx) * nCVI * nHI
              - 10**R_Cm_to_CI_via_HII(Tx) * nCm * nHII
              - 10**R_NI_to_NII_via_HII(Tx) * nNI * nHII
              + 10**R_NII_to_NI_via_HI(Tx) * nNII * nHI
              + 10**R_NIII_to_NII_via_HI(Tx) * nNIII * nHI
              + 10**R_NIV_to_NIII_via_HI(Tx) * nNIV * nHI
              + 10**R_NV_to_NIV_via_HI(Tx) * nNV * nHI
              + 10**R_NVI_to_NV_via_HI(Tx) * nNVI * nHI
              - 10**R_HII_to_HI_via_e_caseA(Tx) * nHII * ne # H CaseA
             )

  dnHm_dt = (
             + 10**R_HI_to_Hm_via_e(Tx) * nHI * ne
             - 10**R_Hm_to_HI_via_HI(Tx) * nHm * nHI
             - 10**R_Hm_to_HI_via_e(Tx) * nHm * ne
             - 10**R_Hm_to_HI_via_HII(Tx) * nHm * nHII
             - 10**R_HeII_to_HeI_via_Hm(Tx) * nHeII * nHm
            )

  dnHeI_dt = (
              - 10**R_HeI_to_HeII_via_HII(Tx) * nHeI * nHII
              + 10**grain_rec_HeII_to_HeI * nHeII * ne # grain_recombination
              - 10**R_HeI_to_HeII_via_e(Tx) * nHeI * ne
              + 10**R_HeII_to_HeI_via_Hm(Tx) * nHeII * nHm
              + 10**R_HeII_to_HeI_via_HI(Tx) * nHeII * nHI
              + 10**R_CI_to_CII_via_HeII(Tx) * nCI * nHeII
              + 10**R_CII_to_CIII_via_HeII(Tx) * nCII * nHeII
              - 10**R_CIV_to_CIII_via_HeI(Tx) * nCIV * nHeI
              - 10**R_CV_to_CIV_via_HeI(Tx) * nCV * nHeI
              + 10**R_NII_to_NIII_via_HeII(Tx) * nNII * nHeII
              - 10**R_NIII_to_NII_via_HeI(Tx) * nNIII * nHeI
              - 10**R_NIV_to_NIII_via_HeI(Tx) * nNIV * nHeI
              - 10**R_NV_to_NIV_via_HeI(Tx) * nNV * nHeI
              + 10**R_HeII_to_HeI_via_e_caseA(Tx) * nHeII * ne # He CaseA
             )

  dnHeII_dt = (
               + 10**R_HeI_to_HeII_via_HII(Tx) * nHeI * nHII
               - 10**grain_rec_HeII_to_HeI * nHeII * ne # grain_recombination
               + 10**R_HeI_to_HeII_via_e(Tx) * nHeI * ne
               - 10**R_HeII_to_HeIII_via_e(Tx) * nHeII * ne
               - 10**R_HeII_to_HeI_via_Hm(Tx) * nHeII * nHm
               - 10**R_HeII_to_HeI_via_HI(Tx) * nHeII * nHI
               + 10**R_HeIII_to_HeII_via_HI(Tx) * nHeIII * nHI
               + 10**R_HeIII_to_HeII_via_e(Tx) * nHeIII * ne
               - 10**R_CI_to_CII_via_HeII(Tx) * nCI * nHeII
               - 10**R_CII_to_CIII_via_HeII(Tx) * nCII * nHeII
               + 10**R_CIV_to_CIII_via_HeI(Tx) * nCIV * nHeI
               + 10**R_CV_to_CIV_via_HeI(Tx) * nCV * nHeI
               - 10**R_NII_to_NIII_via_HeII(Tx) * nNII * nHeII
               + 10**R_NIII_to_NII_via_HeI(Tx) * nNIII * nHeI
               + 10**R_NIV_to_NIII_via_HeI(Tx) * nNIV * nHeI
               + 10**R_NV_to_NIV_via_HeI(Tx) * nNV * nHeI
               - 10**R_HeII_to_HeI_via_e_caseA(Tx) * nHeII * ne # He CaseA
              )

  dnHeIII_dt = (
                + 10**R_HeII_to_HeIII_via_e(Tx) * nHeII * ne
                - 10**R_HeIII_to_HeII_via_HI(Tx) * nHeIII * nHI
                - 10**R_HeIII_to_HeII_via_e(Tx) * nHeIII * ne
               )

  dnCI_dt = (
             - 10**R_CI_to_CII_via_HeII(Tx) * nCI * nHeII
             + 10**grain_rec_CII_to_CI * nCII * ne # grain_recombination
             - 10**R_CI_to_CII_via_HII(Tx) * nCI * nHII
             - 10**R_CI_to_CII_via_e(Tx) * nCI * ne
             + 10**R_CII_to_CI_via_HI(Tx) * nCII * nHI
             + 10**R_CII_to_CI_via_e(Tx) * nCII * ne
             + 10**R_Cm_to_CI_via_HII(Tx) * nCm * nHII
             - const_CI_e_to_Cm_ * nCI * ne # constant rate
            )

  dnCII_dt = (
              + 10**R_CI_to_CII_via_HeII(Tx) * nCI * nHeII
              - 10**grain_rec_CII_to_CI * nCII * ne # grain_recombination
              + 10**R_CI_to_CII_via_HII(Tx) * nCI * nHII
              + 10**R_CI_to_CII_via_e(Tx) * nCI * ne
              - 10**R_CII_to_CI_via_HI(Tx) * nCII * nHI
              - 10**R_CII_to_CIII_via_HeII(Tx) * nCII * nHeII
              - 10**R_CII_to_CI_via_e(Tx) * nCII * ne
              - 10**R_CII_to_CIII_via_e(Tx) * nCII * ne
              + 10**R_CIII_to_CII_via_HI(Tx) * nCIII * nHI
              + 10**R_CIII_to_CII_via_e(Tx) * nCIII * ne
             )

  dnCIII_dt = (
               + 10**R_CII_to_CIII_via_HeII(Tx) * nCII * nHeII
               + 10**R_CII_to_CIII_via_e(Tx) * nCII * ne
               - 10**R_CIII_to_CII_via_HI(Tx) * nCIII * nHI
               - 10**R_CIII_to_CII_via_e(Tx) * nCIII * ne
               - 10**R_CIII_to_CIV_via_e(Tx) * nCIII * ne
               + 10**R_CIV_to_CIII_via_HeI(Tx) * nCIV * nHeI
               + 10**R_CIV_to_CIII_via_HI(Tx) * nCIV * nHI
               + 10**R_CIV_to_CIII_via_e(Tx) * nCIV * ne
              )

  dnCIV_dt = (
              + 10**R_CIII_to_CIV_via_e(Tx) * nCIII * ne
              - 10**R_CIV_to_CIII_via_HeI(Tx) * nCIV * nHeI
              - 10**R_CIV_to_CIII_via_HI(Tx) * nCIV * nHI
              - 10**R_CIV_to_CIII_via_e(Tx) * nCIV * ne
              - 10**R_CIV_to_CV_via_e(Tx) * nCIV * ne
              + 10**R_CV_to_CIV_via_e(Tx) * nCV * ne
              + 10**R_CV_to_CIV_via_HI(Tx) * nCV * nHI
              + 10**R_CV_to_CIV_via_HeI(Tx) * nCV * nHeI
             )

  dnCV_dt = (
             + 10**R_CIV_to_CV_via_e(Tx) * nCIV * ne
             - 10**R_CV_to_CVI_via_e(Tx) * nCV * ne
             - 10**R_CV_to_CIV_via_e(Tx) * nCV * ne
             - 10**R_CV_to_CIV_via_HI(Tx) * nCV * nHI
             - 10**R_CV_to_CIV_via_HeI(Tx) * nCV * nHeI
             + 10**R_CVI_to_CV_via_HI(Tx) * nCVI * nHI
             + 10**R_CVI_to_CV_via_e(Tx) * nCVI * ne
            )

  dnCVI_dt = (
              + 10**R_CV_to_CVI_via_e(Tx) * nCV * ne
              - 10**R_CVI_to_CVII_via_e(Tx) * nCVI * ne
              - 10**R_CVI_to_CV_via_HI(Tx) * nCVI * nHI
              - 10**R_CVI_to_CV_via_e(Tx) * nCVI * ne
              + 10**R_CVII_to_CVI_via_e(Tx) * nCVII * ne
             )

  dnCVII_dt = (
               + 10**R_CVI_to_CVII_via_e(Tx) * nCVI * ne
               - 10**R_CVII_to_CVI_via_e(Tx) * nCVII * ne
              )

  dnCm_dt = (
             - 10**R_Cm_to_CI_via_HII(Tx) * nCm * nHII
             + const_CI_e_to_Cm_ * nCI * ne # constant rate
            )

  dnNI_dt = (
             - 10**R_NI_to_NII_via_HII(Tx) * nNI * nHII
             - 10**R_NI_to_NII_via_e(Tx) * nNI * ne
             + 10**R_NII_to_NI_via_HI(Tx) * nNII * nHI
             + 10**R_NII_to_NI_via_e(Tx) * nNII * ne
            )

  dnNII_dt = (
              + 10**R_NI_to_NII_via_HII(Tx) * nNI * nHII
              + 10**R_NI_to_NII_via_e(Tx) * nNI * ne
              - 10**R_NII_to_NI_via_HI(Tx) * nNII * nHI
              - 10**R_NII_to_NIII_via_HeII(Tx) * nNII * nHeII
              - 10**R_NII_to_NI_via_e(Tx) * nNII * ne
              - 10**R_NII_to_NIII_via_e(Tx) * nNII * ne
              + 10**R_NIII_to_NII_via_HeI(Tx) * nNIII * nHeI
              + 10**R_NIII_to_NII_via_HI(Tx) * nNIII * nHI
              + 10**R_NIII_to_NII_via_e(Tx) * nNIII * ne
             )

  dnNIII_dt = (
               + 10**R_NII_to_NIII_via_HeII(Tx) * nNII * nHeII
               + 10**R_NII_to_NIII_via_e(Tx) * nNII * ne
               - 10**R_NIII_to_NII_via_HeI(Tx) * nNIII * nHeI
               - 10**R_NIII_to_NII_via_HI(Tx) * nNIII * nHI
               - 10**R_NIII_to_NII_via_e(Tx) * nNIII * ne
               - 10**R_NIII_to_NIV_via_e(Tx) * nNIII * ne
               + 10**R_NIV_to_NIII_via_HeI(Tx) * nNIV * nHeI
               + 10**R_NIV_to_NIII_via_HI(Tx) * nNIV * nHI
               + 10**R_NIV_to_NIII_via_e(Tx) * nNIV * ne
              )

  dnNIV_dt = (
              + 10**R_NIII_to_NIV_via_e(Tx) * nNIII * ne
              - 10**R_NIV_to_NIII_via_HeI(Tx) * nNIV * nHeI
              - 10**R_NIV_to_NIII_via_HI(Tx) * nNIV * nHI
              - 10**R_NIV_to_NIII_via_e(Tx) * nNIV * ne
              - 10**R_NIV_to_NV_via_e(Tx) * nNIV * ne
              + 10**R_NV_to_NIV_via_HeI(Tx) * nNV * nHeI
              + 10**R_NV_to_NIV_via_HI(Tx) * nNV * nHI
              + 10**R_NV_to_NIV_via_e(Tx) * nNV * ne
             )

  dnNV_dt = (
             + 10**R_NIV_to_NV_via_e(Tx) * nNIV * ne
             - 10**R_NV_to_NIV_via_HeI(Tx) * nNV * nHeI
             - 10**R_NV_to_NIV_via_HI(Tx) * nNV * nHI
             - 10**R_NV_to_NIV_via_e(Tx) * nNV * ne
             - 10**R_NV_to_NVI_via_e(Tx) * nNV * ne
             + 10**R_NVI_to_NV_via_HI(Tx) * nNVI * nHI
             + 10**R_NVI_to_NV_via_e(Tx) * nNVI * ne
            )

  dnNVI_dt = (
              + 10**R_NV_to_NVI_via_e(Tx) * nNV * ne
              - 10**R_NVI_to_NV_via_HI(Tx) * nNVI * nHI
              - 10**R_NVI_to_NVII_via_e(Tx) * nNVI * ne
              - 10**R_NVI_to_NV_via_e(Tx) * nNVI * ne
              + 10**R_NVII_to_NVI_via_e(Tx) * nNVII * ne
             )

  dnNVII_dt = (
               + 10**R_NVI_to_NVII_via_e(Tx) * nNVI * ne
               - 10**R_NVII_to_NVI_via_e(Tx) * nNVII * ne
               - 10**R_NVII_to_NVIII_via_e(Tx) * nNVII * ne
               + 10**R_NVIII_to_NVII_via_e(Tx) * nNVIII * ne
              )

  dnNVIII_dt = (
                + 10**R_NVII_to_NVIII_via_e(Tx) * nNVII * ne
                - 10**R_NVIII_to_NVII_via_e(Tx) * nNVIII * ne
               )

  Lamb = Lambda(
                T, nHI, nHII, nHm, nHeI, nHeII, nHeIII, nCI, nCII, nCIII, nCIV, 
                nCV, nCVI, nCVII, nCm, nNI, nNII, nNIII, nNIV, nNV, nNVI, 
                nNVII, nNVIII)

  dne_dt = (
           + 1 * dnHII_dt + -1 * dnHm_dt + 1 * dnHeII_dt + 2 * dnHeIII_dt + 1 * dnCII_dt + 2 * dnCIII_dt
           + 3 * dnCIV_dt + 4 * dnCV_dt + 5 * dnCVI_dt + 6 * dnCVII_dt + -1 * dnCm_dt + 1 * dnNII_dt
           + 2 * dnNIII_dt + 3 * dnNIV_dt + 4 * dnNV_dt + 5 * dnNVI_dt + 6 * dnNVII_dt + 7 * dnNVIII_dt
          
           )

  dntot_dt = (
                dne_dt + dnHI_dt + dnHII_dt + dnHm_dt + dnHeI_dt + dnHeII_dt + dnHeIII_dt
              + dnCI_dt + dnCII_dt + dnCIII_dt + dnCIV_dt + dnCV_dt + dnCVI_dt
              + dnCVII_dt + dnCm_dt + dnNI_dt + dnNII_dt + dnNIII_dt + dnNIV_dt
              + dnNV_dt + dnNVI_dt + dnNVII_dt + dnNVIII_dt
             )

  dT_dt = -1.0 * (gamma - 1.0) / kB / ntot * (Lamb + 1. / (gamma - 1.) * kB * T * dntot_dt)

  return [
          dnHI_dt, dnHII_dt, dnHm_dt, dnHeI_dt, dnHeII_dt, dnHeIII_dt,
          dnCI_dt, dnCII_dt, dnCIII_dt, dnCIV_dt, dnCV_dt, dnCVI_dt,
          dnCVII_dt, dnCm_dt, dnNI_dt, dnNII_dt, dnNIII_dt, dnNIV_dt,
          dnNV_dt, dnNVI_dt, dnNVII_dt, dnNVIII_dt, dT_dt
         ]




nH = 1000.0

He_solar = 10**(-1.0) # See Table_1 in Wiersma et al - 2009, 393, 99–107
nHe = He_solar * nH
print('nHe (cm^-3) = ', nHe)

C_solar = 10**(-3.61) # See Table_1 in Wiersma et al - 2009, 393, 99–107
nC = C_solar * nH

N_solar = 10**(-4.07)
nN = N_solar * nH

print('nC (cm^-3) = ', nC)
print()
print('nN (cm^-3) = ', nN)
print()
print('H, C, N before = ', nH, nC, nN)
print()

T_i = 10**6.00

nHm_i = 1e-5 * nH
nH0_i = 0.001 * nH
nHp_i = nH - nH0_i

nHe0_i = 0.0001 * nHe
nHep_i = 0.001 * nHe
nHepp_i= nHe - nHe0_i - nHep_i
 
nCm_i = 1e-6 * nC
nC0_i = 1e-5 * nC
nC1_i = 1e-5 * nC
nC2_i = 1e-4 * nC
nC3_i = 1e-4 * nC
nC4_i = 1e-3 * nC
nC5_i = 1e-2 * nC
nC6_i = nC - (nCm_i + nC0_i + nC1_i + nC2_i + nC3_i + nC4_i + nC5_i)

nN0_i = 1e-5 * nN
nN1_i = 1e-5 * nN
nN2_i = 1e-4 * nN
nN3_i = 1e-4 * nN
nN4_i = 1e-3 * nN
nN5_i = 1e-2 * nN
nN6_i = 1e-2 * nN
nN7_i = nN - (nN0_i + nN1_i + nN2_i + nN3_i + nN4_i + nN5_i + nN6_i)

y0 = [
      nH0_i, nHp_i, nHm_i, nHe0_i, nHep_i, nHepp_i,
      nC0_i, nC1_i, nC2_i, nC3_i, nC4_i, nC5_i, nC6_i, nCm_i,
      nN0_i, nN1_i, nN2_i, nN3_i, nN4_i, nN5_i, nN6_i, nN7_i,
      T_i
      ]


A_v = 1.0
G0 = 0.01
dust_ratio = 1.0

t_span = (1*3.16e7, 3000*3.16e7)

solution = solve_ivp(func, t_span, y0, method='LSODA', dense_output=True)

t = np.linspace(t_span[0], t_span[1], 10000) # This 10000 is not years, it is the number of points in linspace !!!!
y = solution.sol(t)

print(y.shape)

t_yrs = t / 3.16e7

nH0  = y[0, :]
nHp  = y[1, :]
nHm  = y[2, :]
nHe0 = y[3, :]
nHep = y[4, :]
nHepp= y[5, :]

nC0 = y[6, :]
nC1 = y[7, :]
nC2 = y[8, :]
nC3 = y[9, :]
nC4 = y[10, :]
nC5 = y[11, :]
nC6 = y[12, :]
nCm = y[13, :]

nN0 = y[14, :]
nN1 = y[15, :]
nN2 = y[16, :]
nN3 = y[17, :]
nN4 = y[18, :]
nN5 = y[19, :]
nN6 = y[20, :]
nN7 = y[21, :]

T = y[22, :]

print('nC0 = ', nC0)
print()
print('nC1 = ', nC1)
print()
print('nC2 = ', nC2)
print()
print('nC3 = ', nC3)
print()
print('nC4 = ', nC4)
print()
print('nC5 = ', nC5)
print()
print('nC6 = ', nC6)
print()



#------ Result from "test_primordial_hdf5_v2.py" code -----
with open('chimesRes_C_N.pkl', 'rb') as f:
  df = pickle.load(f)
# dictx = {'t_Arr_in_yrs': t_Arr_in_yrs, 'TEvol': TEvol, 'nHe0': nHe0, 'nHep': nHep, 'nHepp': nHepp}
t_Arr_in_yrsx = df['t_Arr_in_yrs']
TEvolx = df['TEvol']
nH0x = df['nH0']
nHpx = df['nHp']
nHx = nH0x + nHpx
nHe0x = df['nHe0']
nHepx = df['nHep']
nHeppx = df['nHepp']
nHeTotx = nHe0x + nHepx + nHeppx

nC0x = df['nC0']
nC1x = df['nC1']
nC2x = df['nC2']
nC3x = df['nC3']
nC4x = df['nC4']
nC5x = df['nC5']
nC6x = df['nC6']
nCx = nC0x + nC1x + nC2x + nC3x + nC4x + nC5x + nC6x

nN0x = df['nN0']
nN1x = df['nN1']
nN2x = df['nN2']
nN3x = df['nN3']
nN4x = df['nN4']
nN5x = df['nN5']
nN6x = df['nN6']
nN7x = df['nN7']
nNx = nN0x + nN1x + nN2x + nN3x + nN4x + nN5x + nN6x + nN7x
#----------------------------------------------------------

ne = nHp + (nHep + 2.0 * nHepp) + (nC1 + 2.0 * nC2 + 3.0 * nC3 + 4.0 * nC4 + 5.0 * nC5 + 6.0 * nC6)
nex = nHpx + (nHepx + 2.0 * nHeppx) + (nC1x + 2.0 * nC2x + 3.0 * nC3x + 4.0 * nC4x + 5.0 * nC5x + 6.0 * nC6x)


plt.figure(figsize = (16, 8))

plt.subplot(2, 3, 1)
plt.scatter(t_yrs, np.log10(T), s = 2, color = 'k', label = 'my own code')
plt.scatter(t_Arr_in_yrsx, np.log10(TEvolx), s = 2, color = 'orange', label = 'chimes result', linestyle = '--')
plt.xlim(0, 3000)
plt.ylim(1, 8)
plt.legend()

plt.subplot(2, 3, 2)

plt.plot(t_yrs, nHe0, color = 'r', label = 'nHe0')
plt.plot(t_yrs, nHep, color = 'g', label = 'nHep')
plt.plot(t_yrs, nHepp, color = 'b', label = 'nHepp')

plt.plot(t_Arr_in_yrsx, nHe0x, color = 'r', label = 'nHe0 - chimes', linestyle = ':')
plt.plot(t_Arr_in_yrsx, nHepx, color = 'g', label = 'nHep - chimes', linestyle = ':')
plt.plot(t_Arr_in_yrsx, nHeppx, color = 'b', label = 'nHepp - chimes', linestyle = ':')

plt.xlim(0, 5000)
plt.ylim(1e-8, 300)

plt.yscale('log')
plt.title('solve_ivp')
plt.legend()


plt.subplot(2, 3, 3)
nHeTot = nHe0 + nHep + nHepp
plt.plot(T, nH0/nH, label = 'nH0', color = 'r')
plt.plot(T, nHp/nH, label = 'nHp', color = 'g')
plt.plot(T, nHm/nH, label = 'nHm', color = 'yellow')
plt.plot(T, nHe0/nHeTot, label = 'nHe0', color = 'b')
plt.plot(T, nHep/nHeTot, label = 'nHep', color = 'orange')
plt.plot(T, nHepp/nHeTot,label = 'nHepp', color = 'purple')

print('nHm/nH = ', np.sort(nHm/nH))

plt.plot(TEvolx, nH0x/nHx, label = 'nH0 - chimes', color = 'r', linestyle = ':')
plt.plot(TEvolx, nHpx/nHx, label = 'nHp - chimes', color = 'g', linestyle = ':')
plt.plot(TEvolx, nHe0x/nHeTotx, label = 'nHe0 - chimes', color = 'b', linestyle = ':')
plt.plot(TEvolx, nHepx/nHeTotx, label = 'nHep - chimes', color = 'orange', linestyle = ':')
plt.plot(TEvolx, nHeppx/nHeTotx,label = 'nHepp - chimes', color = 'purple', linestyle = ':')

plt.yscale('log')
plt.xscale('log')
plt.ylim(2e-10, 1.2)
plt.xlim(1e3, 1e6)
plt.legend()


plt.subplot(2, 3, 4)
plt.plot(T, nC0/nC, label = 'nC0', color = 'r')
plt.plot(T, nC1/nC, label = 'nC1', color = 'g')
plt.plot(T, nC2/nC, label = 'nC2', color = 'b')
plt.plot(T, nC3/nC, label = 'nC3', color = 'orange')
plt.plot(T, nC4/nC, label = 'nC4', color = 'purple')
plt.plot(T, nC5/nC, label = 'nC5', color = 'lime')
plt.plot(T, nC6/nC, label = 'nC6', color = 'pink')

plt.plot(TEvolx, nC0x/nCx, color = 'r', linestyle = ':')
plt.plot(TEvolx, nC1x/nCx, color = 'g', linestyle = ':')
plt.plot(TEvolx, nC2x/nCx, color = 'b', linestyle = ':')
plt.plot(TEvolx, nC3x/nCx, color = 'orange', linestyle = ':')
plt.plot(TEvolx, nC4x/nCx, color = 'purple', linestyle = ':')
plt.plot(TEvolx, nC5x/nCx, color = 'lime', linestyle = ':')
plt.plot(TEvolx, nC6x/nCx, color = 'pink', linestyle = ':')

plt.yscale('log')
plt.xscale('log')
plt.ylim(2e-3, 1.2)
plt.xlim(1e4, 1e6)
#plt.legend()

#***************** Nitrogen *******************
plt.subplot(2, 3, 5)
plt.plot(T, nN0/nN, label = 'nN0', color = 'r')
plt.plot(T, nN1/nN, label = 'nN1', color = 'g')
plt.plot(T, nN2/nN, label = 'nN2', color = 'b')
plt.plot(T, nN3/nN, label = 'nN3', color = 'orange')
plt.plot(T, nN4/nN, label = 'nN4', color = 'purple')
plt.plot(T, nN5/nN, label = 'nN5', color = 'lime')
plt.plot(T, nN6/nN, label = 'nN6', color = 'pink')
#plt.plot(T, nN7/nN, label = 'nN7', color = 'cyan')

plt.plot(TEvolx, nN0x/nNx, label = 'nN0', color = 'r', linestyle = ':')
plt.plot(TEvolx, nN1x/nNx, label = 'nN1', color = 'g', linestyle = ':')
plt.plot(TEvolx, nN2x/nNx, label = 'nN2', color = 'b', linestyle = ':')
plt.plot(TEvolx, nN3x/nNx, label = 'nN3', color = 'orange', linestyle = ':')
plt.plot(TEvolx, nN4x/nNx, label = 'nN4', color = 'purple', linestyle = ':')
plt.plot(TEvolx, nN5x/nNx, label = 'nN5', color = 'lime', linestyle = ':')
plt.plot(TEvolx, nN6x/nNx, label = 'nN6', color = 'pink', linestyle = ':')
plt.plot(TEvolx, nN7x/nNx, label = 'nN7', color = 'cyan', linestyle = ':', linewidth = 5)

plt.yscale('log')
plt.xscale('log')
plt.ylim(2e-3, 1.2)
plt.xlim(1e4, 1e6)
#plt.legend()

plt.tight_layout()

plt.savefig('result_v5_HCN.png')

plt.show()



