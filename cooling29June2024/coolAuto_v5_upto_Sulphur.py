
import numpy as np
import matplotlib.pyplot as plt
#from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp
import pickle
from datax_vz import *
from scipy.interpolate import RegularGridInterpolator
import time


#----- Lambda
def Lambda(T, nHI, nHII, nHm, nHeI, nHeII, nHeIII, nCI, nCII, nCIII, nCIV, nCV, 
           nCVI, nCVII, nCm, nNI, nNII, nNIII, nNIV, nNV, nNVI, nNVII, 
           nNVIII, nOI, nOII, nOIII, nOIV, nOV, nOVI, nOVII, nOVIII, nOIX, 
           nOm, nNeI, nNeII, nNeIII, nNeIV, nNeV, nNeVI, nNeVII, nNeVIII, nNeIX, 
           nNeX, nNeXI, nMgI, nMgII, nMgIII, nMgIV, nMgV, nMgVI, nMgVII, nMgVIII, 
           nMgIX, nMgX, nMgXI, nMgXII, nMgXIII, nSiI, nSiII, nSiIII, nSiIV, nSiV, 
           nSiVI, nSiVII, nSiVIII, nSiIX, nSiX, nSiXI, nSiXII, nSiXIII, nSiXIV, nSiXV, 
           nSI, nSII, nSIII, nSIV, nSV, nSVI, nSVII, nSVIII, nSIX, nSX, 
           nSXI, nSXII, nSXIII, nSXIV, nSXV, nSXVI, nSXVII):

  Tx = np.log10(T)

  ne = (
         1 * nHII - nHm + (nHeII + 2.0 * nHeIII) + 1 * nCII + 2 * nCIII + 3 * nCIV
       + 4 * nCV + 5 * nCVI + 6 * nCVII + -1 * nCm + 1 * nNII
       + 2 * nNIII + 3 * nNIV + 4 * nNV + 5 * nNVI + 6 * nNVII + 7 * nNVIII
       + 1 * nOII + 2 * nOIII + 3 * nOIV + 4 * nOV + 5 * nOVI
       + 6 * nOVII + 7 * nOVIII + 8 * nOIX + -1 * nOm + 1 * nNeII
       + 2 * nNeIII + 3 * nNeIV + 4 * nNeV + 5 * nNeVI + 6 * nNeVII + 7 * nNeVIII
       + 8 * nNeIX + 9 * nNeX + 10 * nNeXI + 1 * nMgII + 2 * nMgIII
       + 3 * nMgIV + 4 * nMgV + 5 * nMgVI + 6 * nMgVII + 7 * nMgVIII + 8 * nMgIX
       + 9 * nMgX + 10 * nMgXI + 11 * nMgXII + 12 * nMgXIII + 1 * nSiII
       + 2 * nSiIII + 3 * nSiIV + 4 * nSiV + 5 * nSiVI + 6 * nSiVII + 7 * nSiVIII
       + 8 * nSiIX + 9 * nSiX + 10 * nSiXI + 11 * nSiXII + 12 * nSiXIII + 13 * nSiXIV
       + 14 * nSiXV + 1 * nSII + 2 * nSIII + 3 * nSIV + 4 * nSV
       + 5 * nSVI + 6 * nSVII + 7 * nSVIII + 8 * nSIX + 9 * nSX + 10 * nSXI
       + 11 * nSXII + 12 * nSXIII + 13 * nSXIV + 14 * nSXV + 15 * nSXVI + 16 * nSXVII
      
       )

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
        + 10**cooling_rate_4d("OI", T, nHI, ne, nHII, Temp_4d, HIDensity_4d, elecDensity_4d, HIIDensity_4d) * nOI * ne # cooling via OI
        + 10**cooling_rate_2d("SiII", T, ne, Temp_2d, elecDensity_2d) * nSiII * ne # cooling via SiII
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
        + 10**gOII(Tx) * nOII * ne
        + 10**grain_recomb_cooling_rate * nOII * ne # grain_recombination cooling!
        + 10**gOIII(Tx) * nOIII * ne
        + 10**gOIV(Tx) * nOIV * ne
        + 10**gOV(Tx) * nOV * ne
        + 10**gOVI(Tx) * nOVI * ne
        + 10**gOVII(Tx) * nOVII * ne
        + 10**gOVIII(Tx) * nOVIII * ne
        + 10**gOIX(Tx) * nOIX * ne
        + 10**gNeI(Tx) * nNeI * ne
        + 10**gNeII(Tx) * nNeII * ne
        + 10**gNeIII(Tx) * nNeIII * ne
        + 10**gNeIV(Tx) * nNeIV * ne
        + 10**gNeV(Tx) * nNeV * ne
        + 10**gNeVI(Tx) * nNeVI * ne
        + 10**gNeVII(Tx) * nNeVII * ne
        + 10**gNeVIII(Tx) * nNeVIII * ne
        + 10**gNeIX(Tx) * nNeIX * ne
        + 10**gNeX(Tx) * nNeX * ne
        + 10**gNeXI(Tx) * nNeXI * ne
        + 10**gMgI(Tx) * nMgI * ne
        + 10**gMgII(Tx) * nMgII * ne
        + 10**grain_recomb_cooling_rate * nMgII * ne # grain_recombination cooling!
        + 10**gMgIII(Tx) * nMgIII * ne
        + 10**gMgIV(Tx) * nMgIV * ne
        + 10**gMgV(Tx) * nMgV * ne
        + 10**gMgVI(Tx) * nMgVI * ne
        + 10**gMgVII(Tx) * nMgVII * ne
        + 10**gMgVIII(Tx) * nMgVIII * ne
        + 10**gMgIX(Tx) * nMgIX * ne
        + 10**gMgX(Tx) * nMgX * ne
        + 10**gMgXI(Tx) * nMgXI * ne
        + 10**gMgXII(Tx) * nMgXII * ne
        + 10**gMgXIII(Tx) * nMgXIII * ne
        + 10**gSiI(Tx) * nSiI * ne
        + 10**grain_recomb_cooling_rate * nSiII * ne # grain_recombination cooling!
        + 10**gSiIII(Tx) * nSiIII * ne
        + 10**gSiIV(Tx) * nSiIV * ne
        + 10**gSiV(Tx) * nSiV * ne
        + 10**gSiVI(Tx) * nSiVI * ne
        + 10**gSiVII(Tx) * nSiVII * ne
        + 10**gSiVIII(Tx) * nSiVIII * ne
        + 10**gSiIX(Tx) * nSiIX * ne
        + 10**gSiX(Tx) * nSiX * ne
        + 10**gSiXI(Tx) * nSiXI * ne
        + 10**gSiXII(Tx) * nSiXII * ne
        + 10**gSiXIII(Tx) * nSiXIII * ne
        + 10**gSiXIV(Tx) * nSiXIV * ne
        + 10**gSiXV(Tx) * nSiXV * ne
        + 10**gSI(Tx) * nSI * ne
        + 10**gSII(Tx) * nSII * ne
        + 10**grain_recomb_cooling_rate * nSII * ne # grain_recombination cooling!
        + 10**gSIII(Tx) * nSIII * ne
        + 10**gSIV(Tx) * nSIV * ne
        + 10**gSV(Tx) * nSV * ne
        + 10**gSVI(Tx) * nSVI * ne
        + 10**gSVII(Tx) * nSVII * ne
        + 10**gSVIII(Tx) * nSVIII * ne
        + 10**gSIX(Tx) * nSIX * ne
        + 10**gSX(Tx) * nSX * ne
        + 10**gSXI(Tx) * nSXI * ne
        + 10**gSXII(Tx) * nSXII * ne
        + 10**gSXIII(Tx) * nSXIII * ne
        + 10**gSXIV(Tx) * nSXIV * ne
        + 10**gSXV(Tx) * nSXV * ne
        + 10**gSXVI(Tx) * nSXVI * ne
        + 10**gSXVII(Tx) * nSXVII * ne
        + LCompton)

  return Lamb



#----- func
def func(t, y):

  nHI, nHII, nHm, nHeI, nHeII, nHeIII, nCI, nCII, nCIII, nCIV, \
  nCV, nCVI, nCVII, nCm, nNI, nNII, nNIII, nNIV, nNV, nNVI, \
  nNVII, nNVIII, nOI, nOII, nOIII, nOIV, nOV, nOVI, nOVII, nOVIII, \
  nOIX, nOm, nNeI, nNeII, nNeIII, nNeIV, nNeV, nNeVI, nNeVII, nNeVIII, \
  nNeIX, nNeX, nNeXI, nMgI, nMgII, nMgIII, nMgIV, nMgV, nMgVI, nMgVII, \
  nMgVIII, nMgIX, nMgX, nMgXI, nMgXII, nMgXIII, nSiI, nSiII, nSiIII, nSiIV, \
  nSiV, nSiVI, nSiVII, nSiVIII, nSiIX, nSiX, nSiXI, nSiXII, nSiXIII, nSiXIV, \
  nSiXV, nSI, nSII, nSIII, nSIV, nSV, nSVI, nSVII, nSVIII, nSIX, \
  nSX, nSXI, nSXII, nSXIII, nSXIV, nSXV, nSXVI, nSXVII, T = y

  Tx = np.log10(T)

  ne = (
       + 1 * nHII + -1 * nHm + 1 * nHeII + 2 * nHeIII + 1 * nCII + 2 * nCIII
       + 3 * nCIV + 4 * nCV + 5 * nCVI + 6 * nCVII + -1 * nCm + 1 * nNII
       + 2 * nNIII + 3 * nNIV + 4 * nNV + 5 * nNVI + 6 * nNVII + 7 * nNVIII
       + 1 * nOII + 2 * nOIII + 3 * nOIV + 4 * nOV + 5 * nOVI + 6 * nOVII
       + 7 * nOVIII + 8 * nOIX + -1 * nOm + 1 * nNeII + 2 * nNeIII + 3 * nNeIV
       + 4 * nNeV + 5 * nNeVI + 6 * nNeVII + 7 * nNeVIII + 8 * nNeIX + 9 * nNeX
       + 10 * nNeXI + 1 * nMgII + 2 * nMgIII + 3 * nMgIV + 4 * nMgV + 5 * nMgVI
       + 6 * nMgVII + 7 * nMgVIII + 8 * nMgIX + 9 * nMgX + 10 * nMgXI + 11 * nMgXII
       + 12 * nMgXIII + 1 * nSiII + 2 * nSiIII + 3 * nSiIV + 4 * nSiV + 5 * nSiVI
       + 6 * nSiVII + 7 * nSiVIII + 8 * nSiIX + 9 * nSiX + 10 * nSiXI + 11 * nSiXII
       + 12 * nSiXIII + 13 * nSiXIV + 14 * nSiXV + 1 * nSII + 2 * nSIII + 3 * nSIV
       + 4 * nSV + 5 * nSVI + 6 * nSVII + 7 * nSVIII + 8 * nSIX + 9 * nSX
       + 10 * nSXI + 11 * nSXII + 12 * nSXIII + 13 * nSXIV + 14 * nSXV + 15 * nSXVI
       + 16 * nSXVII
       )

  ntot = (
           ne + nHI + nHII + nHm + nHeI + nHeII + nHeIII
         + nCI + nCII + nCIII + nCIV + nCV + nCVI
         + nCVII + nCm + nNI + nNII + nNIII + nNIV
         + nNV + nNVI + nNVII + nNVIII + nOI + nOII
         + nOIII + nOIV + nOV + nOVI + nOVII + nOVIII
         + nOIX + nOm + nNeI + nNeII + nNeIII + nNeIV
         + nNeV + nNeVI + nNeVII + nNeVIII + nNeIX + nNeX
         + nNeXI + nMgI + nMgII + nMgIII + nMgIV + nMgV
         + nMgVI + nMgVII + nMgVIII + nMgIX + nMgX + nMgXI
         + nMgXII + nMgXIII + nSiI + nSiII + nSiIII + nSiIV
         + nSiV + nSiVI + nSiVII + nSiVIII + nSiIX + nSiX
         + nSiXI + nSiXII + nSiXIII + nSiXIV + nSiXV + nSI
         + nSII + nSIII + nSIV + nSV + nSVI + nSVII
         + nSVIII + nSIX + nSX + nSXI + nSXII + nSXIII
         + nSXIV + nSXV + nSXVI + nSXVII
         )

  grain_rec_HII_to_HI = grain_recomb_rate("HII", T, ne, G0, A_v, Temp, Psi)
  grain_rec_HeII_to_HeI = grain_recomb_rate("HeII", T, ne, G0, A_v, Temp, Psi)
  grain_rec_CII_to_CI = grain_recomb_rate("CII", T, ne, G0, A_v, Temp, Psi)
  grain_rec_OII_to_OI = grain_recomb_rate("OII", T, ne, G0, A_v, Temp, Psi)
  grain_rec_SiII_to_SiI = grain_recomb_rate("SiII", T, ne, G0, A_v, Temp, Psi)
  grain_rec_MgII_to_MgI = grain_recomb_rate("MgII", T, ne, G0, A_v, Temp, Psi)
  grain_rec_SII_to_SI = grain_recomb_rate("SII", T, ne, G0, A_v, Temp, Psi)

  const_CII_SiI_to_CI_SiII = 2.1000E-09 # constant/rates
  const_OI_e_to_Om_ = 1.5000E-15 # constant/rates
  const_CII_MgI_to_CI_MgII = 1.1000E-09 # constant/rates
  const_NII_MgI_to_NI_MgII = 1.2000E-09 # constant/rates
  const_CI_e_to_Cm_ = 2.2500E-15 # constant/rates
  const_SII_MgI_to_SI_MgII = 2.8000E-10 # constant/rates
  const_SiII_MgI_to_SiI_MgII = 2.9000E-09 # constant/rates

  dnHI_dt = (
             - 10**R_HI_to_HII_via_e(Tx) * nHI * ne
             + 10**grain_rec_HII_to_HI * nHII * ne # grain_recombination
             - 10**R_HI_to_Hm_via_e(Tx) * nHI * ne
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
             + 10**R_OI_to_OII_via_HII(Tx) * nOI * nHII
             - 10**R_OII_to_OI_via_HI(Tx) * nOII * nHI
             - 10**R_OIII_to_OII_via_HI(Tx) * nOIII * nHI
             - 10**R_OIV_to_OIII_via_HI(Tx) * nOIV * nHI
             - 10**R_OV_to_OIV_via_HI(Tx) * nOV * nHI
             - 10**R_OVI_to_OV_via_HI(Tx) * nOVI * nHI
             + 10**R_Om_to_OI_via_HII(Tx) * nOm * nHII
             - 10**R_NeIII_to_NeII_via_HI(Tx) * nNeIII * nHI
             - 10**R_NeIV_to_NeIII_via_HI(Tx) * nNeIV * nHI
             - 10**R_NeV_to_NeIV_via_HI(Tx) * nNeV * nHI
             - 10**R_NeVI_to_NeV_via_HI(Tx) * nNeVI * nHI
             + 10**R_MgI_to_MgII_via_HII(Tx) * nMgI * nHII
             + 10**R_MgII_to_MgIII_via_HII(Tx) * nMgII * nHII
             - 10**R_MgIII_to_MgII_via_HI(Tx) * nMgIII * nHI
             - 10**R_MgIV_to_MgIII_via_HI(Tx) * nMgIV * nHI
             - 10**R_MgV_to_MgIV_via_HI(Tx) * nMgV * nHI
             - 10**R_MgVI_to_MgV_via_HI(Tx) * nMgVI * nHI
             + 10**R_SiI_to_SiII_via_HII(Tx) * nSiI * nHII
             + 10**R_SiII_to_SiIII_via_HII(Tx) * nSiII * nHII
             - 10**R_SiIII_to_SiII_via_HI(Tx) * nSiIII * nHI
             - 10**R_SiIV_to_SiIII_via_HI(Tx) * nSiIV * nHI
             - 10**R_SiV_to_SiIV_via_HI(Tx) * nSiV * nHI
             - 10**R_SiVI_to_SiV_via_HI(Tx) * nSiVI * nHI
             + 10**R_SI_to_SII_via_HII(Tx) * nSI * nHII
             - 10**R_SII_to_SI_via_HI(Tx) * nSII * nHI
             - 10**R_SIII_to_SII_via_HI(Tx) * nSIII * nHI
             - 10**R_SIV_to_SIII_via_HI(Tx) * nSIV * nHI
             - 10**R_SV_to_SIV_via_HI(Tx) * nSV * nHI
             - 10**R_SVI_to_SV_via_HI(Tx) * nSVI * nHI
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
              - 10**R_OI_to_OII_via_HII(Tx) * nOI * nHII
              + 10**R_OII_to_OI_via_HI(Tx) * nOII * nHI
              + 10**R_OIII_to_OII_via_HI(Tx) * nOIII * nHI
              + 10**R_OIV_to_OIII_via_HI(Tx) * nOIV * nHI
              + 10**R_OV_to_OIV_via_HI(Tx) * nOV * nHI
              + 10**R_OVI_to_OV_via_HI(Tx) * nOVI * nHI
              - 10**R_Om_to_OI_via_HII(Tx) * nOm * nHII
              + 10**R_NeIII_to_NeII_via_HI(Tx) * nNeIII * nHI
              + 10**R_NeIV_to_NeIII_via_HI(Tx) * nNeIV * nHI
              + 10**R_NeV_to_NeIV_via_HI(Tx) * nNeV * nHI
              + 10**R_NeVI_to_NeV_via_HI(Tx) * nNeVI * nHI
              - 10**R_MgI_to_MgII_via_HII(Tx) * nMgI * nHII
              - 10**R_MgII_to_MgIII_via_HII(Tx) * nMgII * nHII
              + 10**R_MgIII_to_MgII_via_HI(Tx) * nMgIII * nHI
              + 10**R_MgIV_to_MgIII_via_HI(Tx) * nMgIV * nHI
              + 10**R_MgV_to_MgIV_via_HI(Tx) * nMgV * nHI
              + 10**R_MgVI_to_MgV_via_HI(Tx) * nMgVI * nHI
              - 10**R_SiI_to_SiII_via_HII(Tx) * nSiI * nHII
              - 10**R_SiII_to_SiIII_via_HII(Tx) * nSiII * nHII
              + 10**R_SiIII_to_SiII_via_HI(Tx) * nSiIII * nHI
              + 10**R_SiIV_to_SiIII_via_HI(Tx) * nSiIV * nHI
              + 10**R_SiV_to_SiIV_via_HI(Tx) * nSiV * nHI
              + 10**R_SiVI_to_SiV_via_HI(Tx) * nSiVI * nHI
              - 10**R_SI_to_SII_via_HII(Tx) * nSI * nHII
              + 10**R_SII_to_SI_via_HI(Tx) * nSII * nHI
              + 10**R_SIII_to_SII_via_HI(Tx) * nSIII * nHI
              + 10**R_SIV_to_SIII_via_HI(Tx) * nSIV * nHI
              + 10**R_SV_to_SIV_via_HI(Tx) * nSV * nHI
              + 10**R_SVI_to_SV_via_HI(Tx) * nSVI * nHI
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
              + 10**R_OI_to_OII_via_HeII(Tx) * nOI * nHeII
              - 10**R_OIII_to_OII_via_HeI(Tx) * nOIII * nHeI
              - 10**R_OIV_to_OIII_via_HeI(Tx) * nOIV * nHeI
              - 10**R_OV_to_OIV_via_HeI(Tx) * nOV * nHeI
              - 10**R_NeIII_to_NeII_via_HeI(Tx) * nNeIII * nHeI
              - 10**R_NeIV_to_NeIII_via_HeI(Tx) * nNeIV * nHeI
              - 10**R_NeV_to_NeIV_via_HeI(Tx) * nNeV * nHeI
              - 10**R_MgIV_to_MgIII_via_HeI(Tx) * nMgIV * nHeI
              - 10**R_MgV_to_MgIV_via_HeI(Tx) * nMgV * nHeI
              + 10**R_SiI_to_SiII_via_HeII(Tx) * nSiI * nHeII
              + 10**R_SiII_to_SiIII_via_HeII(Tx) * nSiII * nHeII
              + 10**R_SiIII_to_SiIV_via_HeII(Tx) * nSiIII * nHeII
              - 10**R_SiIV_to_SiIII_via_HeI(Tx) * nSiIV * nHeI
              - 10**R_SiV_to_SiIV_via_HeI(Tx) * nSiV * nHeI
              + 10**R_SII_to_SIII_via_HeII(Tx) * nSII * nHeII
              + 10**R_SIII_to_SIV_via_HeII(Tx) * nSIII * nHeII
              - 10**R_SIV_to_SIII_via_HeI(Tx) * nSIV * nHeI
              - 10**R_SV_to_SIV_via_HeI(Tx) * nSV * nHeI
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
               - 10**R_OI_to_OII_via_HeII(Tx) * nOI * nHeII
               + 10**R_OIII_to_OII_via_HeI(Tx) * nOIII * nHeI
               + 10**R_OIV_to_OIII_via_HeI(Tx) * nOIV * nHeI
               + 10**R_OV_to_OIV_via_HeI(Tx) * nOV * nHeI
               + 10**R_NeIII_to_NeII_via_HeI(Tx) * nNeIII * nHeI
               + 10**R_NeIV_to_NeIII_via_HeI(Tx) * nNeIV * nHeI
               + 10**R_NeV_to_NeIV_via_HeI(Tx) * nNeV * nHeI
               + 10**R_MgIV_to_MgIII_via_HeI(Tx) * nMgIV * nHeI
               + 10**R_MgV_to_MgIV_via_HeI(Tx) * nMgV * nHeI
               - 10**R_SiI_to_SiII_via_HeII(Tx) * nSiI * nHeII
               - 10**R_SiII_to_SiIII_via_HeII(Tx) * nSiII * nHeII
               - 10**R_SiIII_to_SiIV_via_HeII(Tx) * nSiIII * nHeII
               + 10**R_SiIV_to_SiIII_via_HeI(Tx) * nSiIV * nHeI
               + 10**R_SiV_to_SiIV_via_HeI(Tx) * nSiV * nHeI
               - 10**R_SII_to_SIII_via_HeII(Tx) * nSII * nHeII
               - 10**R_SIII_to_SIV_via_HeII(Tx) * nSIII * nHeII
               + 10**R_SIV_to_SIII_via_HeI(Tx) * nSIV * nHeI
               + 10**R_SV_to_SIV_via_HeI(Tx) * nSV * nHeI
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
             + const_CII_SiI_to_CI_SiII * nCII * nSiI # constant rate
             + const_CII_MgI_to_CI_MgII * nCII * nMgI # constant rate
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
              - const_CII_SiI_to_CI_SiII * nCII * nSiI # constant rate
              - const_CII_MgI_to_CI_MgII * nCII * nMgI # constant rate
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
             + const_NII_MgI_to_NI_MgII * nNII * nMgI # constant rate
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
              - const_NII_MgI_to_NI_MgII * nNII * nMgI # constant rate
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

  dnOI_dt = (
             - 10**R_OI_to_OII_via_HeII(Tx) * nOI * nHeII
             + 10**grain_rec_OII_to_OI * nOII * ne # grain_recombination
             - 10**R_OI_to_OII_via_e(Tx) * nOI * ne
             - 10**R_OI_to_OII_via_HII(Tx) * nOI * nHII
             + 10**R_OII_to_OI_via_HI(Tx) * nOII * nHI
             + 10**R_OII_to_OI_via_e(Tx) * nOII * ne
             + 10**R_Om_to_OI_via_HII(Tx) * nOm * nHII
             - const_OI_e_to_Om_ * nOI * ne # constant rate
            )

  dnOII_dt = (
              + 10**R_OI_to_OII_via_HeII(Tx) * nOI * nHeII
              - 10**grain_rec_OII_to_OI * nOII * ne # grain_recombination
              + 10**R_OI_to_OII_via_e(Tx) * nOI * ne
              + 10**R_OI_to_OII_via_HII(Tx) * nOI * nHII
              - 10**R_OII_to_OI_via_HI(Tx) * nOII * nHI
              - 10**R_OII_to_OI_via_e(Tx) * nOII * ne
              - 10**R_OII_to_OIII_via_e(Tx) * nOII * ne
              + 10**R_OIII_to_OII_via_HeI(Tx) * nOIII * nHeI
              + 10**R_OIII_to_OII_via_HI(Tx) * nOIII * nHI
              + 10**R_OIII_to_OII_via_e(Tx) * nOIII * ne
             )

  dnOIII_dt = (
               + 10**R_OII_to_OIII_via_e(Tx) * nOII * ne
               - 10**R_OIII_to_OII_via_HeI(Tx) * nOIII * nHeI
               - 10**R_OIII_to_OII_via_HI(Tx) * nOIII * nHI
               - 10**R_OIII_to_OII_via_e(Tx) * nOIII * ne
               - 10**R_OIII_to_OIV_via_e(Tx) * nOIII * ne
               + 10**R_OIV_to_OIII_via_e(Tx) * nOIV * ne
               + 10**R_OIV_to_OIII_via_HI(Tx) * nOIV * nHI
               + 10**R_OIV_to_OIII_via_HeI(Tx) * nOIV * nHeI
              )

  dnOIV_dt = (
              + 10**R_OIII_to_OIV_via_e(Tx) * nOIII * ne
              - 10**R_OIV_to_OV_via_e(Tx) * nOIV * ne
              - 10**R_OIV_to_OIII_via_e(Tx) * nOIV * ne
              - 10**R_OIV_to_OIII_via_HI(Tx) * nOIV * nHI
              - 10**R_OIV_to_OIII_via_HeI(Tx) * nOIV * nHeI
              + 10**R_OV_to_OIV_via_HeI(Tx) * nOV * nHeI
              + 10**R_OV_to_OIV_via_HI(Tx) * nOV * nHI
              + 10**R_OV_to_OIV_via_e(Tx) * nOV * ne
             )

  dnOV_dt = (
             + 10**R_OIV_to_OV_via_e(Tx) * nOIV * ne
             - 10**R_OV_to_OIV_via_HeI(Tx) * nOV * nHeI
             - 10**R_OV_to_OIV_via_HI(Tx) * nOV * nHI
             - 10**R_OV_to_OIV_via_e(Tx) * nOV * ne
             - 10**R_OV_to_OVI_via_e(Tx) * nOV * ne
             + 10**R_OVI_to_OV_via_HI(Tx) * nOVI * nHI
             + 10**R_OVI_to_OV_via_e(Tx) * nOVI * ne
            )

  dnOVI_dt = (
              + 10**R_OV_to_OVI_via_e(Tx) * nOV * ne
              - 10**R_OVI_to_OV_via_HI(Tx) * nOVI * nHI
              - 10**R_OVI_to_OV_via_e(Tx) * nOVI * ne
              - 10**R_OVI_to_OVII_via_e(Tx) * nOVI * ne
              + 10**R_OVII_to_OVI_via_e(Tx) * nOVII * ne
             )

  dnOVII_dt = (
               + 10**R_OVI_to_OVII_via_e(Tx) * nOVI * ne
               - 10**R_OVII_to_OVI_via_e(Tx) * nOVII * ne
               - 10**R_OVII_to_OVIII_via_e(Tx) * nOVII * ne
               + 10**R_OVIII_to_OVII_via_e(Tx) * nOVIII * ne
              )

  dnOVIII_dt = (
                + 10**R_OVII_to_OVIII_via_e(Tx) * nOVII * ne
                - 10**R_OVIII_to_OVII_via_e(Tx) * nOVIII * ne
                - 10**R_OVIII_to_OIX_via_e(Tx) * nOVIII * ne
                + 10**R_OIX_to_OVIII_via_e(Tx) * nOIX * ne
               )

  dnOIX_dt = (
              + 10**R_OVIII_to_OIX_via_e(Tx) * nOVIII * ne
              - 10**R_OIX_to_OVIII_via_e(Tx) * nOIX * ne
             )

  dnOm_dt = (
             - 10**R_Om_to_OI_via_HII(Tx) * nOm * nHII
             + const_OI_e_to_Om_ * nOI * ne # constant rate
            )

  dnNeI_dt = (
              - 10**R_NeI_to_NeII_via_e(Tx) * nNeI * ne
              + 10**R_NeII_to_NeI_via_e(Tx) * nNeII * ne
             )

  dnNeII_dt = (
               + 10**R_NeI_to_NeII_via_e(Tx) * nNeI * ne
               - 10**R_NeII_to_NeI_via_e(Tx) * nNeII * ne
               - 10**R_NeII_to_NeIII_via_e(Tx) * nNeII * ne
               + 10**R_NeIII_to_NeII_via_e(Tx) * nNeIII * ne
               + 10**R_NeIII_to_NeII_via_HI(Tx) * nNeIII * nHI
               + 10**R_NeIII_to_NeII_via_HeI(Tx) * nNeIII * nHeI
              )

  dnNeIII_dt = (
                + 10**R_NeII_to_NeIII_via_e(Tx) * nNeII * ne
                - 10**R_NeIII_to_NeII_via_e(Tx) * nNeIII * ne
                - 10**R_NeIII_to_NeIV_via_e(Tx) * nNeIII * ne
                - 10**R_NeIII_to_NeII_via_HI(Tx) * nNeIII * nHI
                - 10**R_NeIII_to_NeII_via_HeI(Tx) * nNeIII * nHeI
                + 10**R_NeIV_to_NeIII_via_e(Tx) * nNeIV * ne
                + 10**R_NeIV_to_NeIII_via_HeI(Tx) * nNeIV * nHeI
                + 10**R_NeIV_to_NeIII_via_HI(Tx) * nNeIV * nHI
               )

  dnNeIV_dt = (
               + 10**R_NeIII_to_NeIV_via_e(Tx) * nNeIII * ne
               - 10**R_NeIV_to_NeIII_via_e(Tx) * nNeIV * ne
               - 10**R_NeIV_to_NeV_via_e(Tx) * nNeIV * ne
               - 10**R_NeIV_to_NeIII_via_HeI(Tx) * nNeIV * nHeI
               - 10**R_NeIV_to_NeIII_via_HI(Tx) * nNeIV * nHI
               + 10**R_NeV_to_NeIV_via_HeI(Tx) * nNeV * nHeI
               + 10**R_NeV_to_NeIV_via_HI(Tx) * nNeV * nHI
               + 10**R_NeV_to_NeIV_via_e(Tx) * nNeV * ne
              )

  dnNeV_dt = (
              + 10**R_NeIV_to_NeV_via_e(Tx) * nNeIV * ne
              - 10**R_NeV_to_NeIV_via_HeI(Tx) * nNeV * nHeI
              - 10**R_NeV_to_NeIV_via_HI(Tx) * nNeV * nHI
              - 10**R_NeV_to_NeIV_via_e(Tx) * nNeV * ne
              - 10**R_NeV_to_NeVI_via_e(Tx) * nNeV * ne
              + 10**R_NeVI_to_NeV_via_HI(Tx) * nNeVI * nHI
              + 10**R_NeVI_to_NeV_via_e(Tx) * nNeVI * ne
             )

  dnNeVI_dt = (
               + 10**R_NeV_to_NeVI_via_e(Tx) * nNeV * ne
               - 10**R_NeVI_to_NeV_via_HI(Tx) * nNeVI * nHI
               - 10**R_NeVI_to_NeV_via_e(Tx) * nNeVI * ne
               - 10**R_NeVI_to_NeVII_via_e(Tx) * nNeVI * ne
               + 10**R_NeVII_to_NeVI_via_e(Tx) * nNeVII * ne
              )

  dnNeVII_dt = (
                + 10**R_NeVI_to_NeVII_via_e(Tx) * nNeVI * ne
                - 10**R_NeVII_to_NeVI_via_e(Tx) * nNeVII * ne
                - 10**R_NeVII_to_NeVIII_via_e(Tx) * nNeVII * ne
                + 10**R_NeVIII_to_NeVII_via_e(Tx) * nNeVIII * ne
               )

  dnNeVIII_dt = (
                 + 10**R_NeVII_to_NeVIII_via_e(Tx) * nNeVII * ne
                 - 10**R_NeVIII_to_NeVII_via_e(Tx) * nNeVIII * ne
                 - 10**R_NeVIII_to_NeIX_via_e(Tx) * nNeVIII * ne
                 + 10**R_NeIX_to_NeVIII_via_e(Tx) * nNeIX * ne
                )

  dnNeIX_dt = (
               + 10**R_NeVIII_to_NeIX_via_e(Tx) * nNeVIII * ne
               - 10**R_NeIX_to_NeVIII_via_e(Tx) * nNeIX * ne
               - 10**R_NeIX_to_NeX_via_e(Tx) * nNeIX * ne
               + 10**R_NeX_to_NeIX_via_e(Tx) * nNeX * ne
              )

  dnNeX_dt = (
              + 10**R_NeIX_to_NeX_via_e(Tx) * nNeIX * ne
              - 10**R_NeX_to_NeXI_via_e(Tx) * nNeX * ne
              - 10**R_NeX_to_NeIX_via_e(Tx) * nNeX * ne
              + 10**R_NeXI_to_NeX_via_e(Tx) * nNeXI * ne
             )

  dnNeXI_dt = (
               + 10**R_NeX_to_NeXI_via_e(Tx) * nNeX * ne
               - 10**R_NeXI_to_NeX_via_e(Tx) * nNeXI * ne
              )

  dnMgI_dt = (
              - 10**R_MgI_to_MgII_via_HII(Tx) * nMgI * nHII
              + 10**grain_rec_MgII_to_MgI * nMgII * ne # grain_recombination
              - 10**R_MgI_to_MgII_via_e(Tx) * nMgI * ne
              + 10**R_MgII_to_MgI_via_e(Tx) * nMgII * ne
              - const_CII_MgI_to_CI_MgII * nCII * nMgI # constant rate
              - const_NII_MgI_to_NI_MgII * nNII * nMgI # constant rate
              - const_SII_MgI_to_SI_MgII * nSII * nMgI # constant rate
              - const_SiII_MgI_to_SiI_MgII * nSiII * nMgI # constant rate
             )

  dnMgII_dt = (
               + 10**R_MgI_to_MgII_via_HII(Tx) * nMgI * nHII
               - 10**grain_rec_MgII_to_MgI * nMgII * ne # grain_recombination
               + 10**R_MgI_to_MgII_via_e(Tx) * nMgI * ne
               - 10**R_MgII_to_MgIII_via_HII(Tx) * nMgII * nHII
               - 10**R_MgII_to_MgI_via_e(Tx) * nMgII * ne
               - 10**R_MgII_to_MgIII_via_e(Tx) * nMgII * ne
               + 10**R_MgIII_to_MgII_via_HI(Tx) * nMgIII * nHI
               + 10**R_MgIII_to_MgII_via_e(Tx) * nMgIII * ne
               + const_CII_MgI_to_CI_MgII * nCII * nMgI # constant rate
               + const_NII_MgI_to_NI_MgII * nNII * nMgI # constant rate
               + const_SII_MgI_to_SI_MgII * nSII * nMgI # constant rate
               + const_SiII_MgI_to_SiI_MgII * nSiII * nMgI # constant rate
              )

  dnMgIII_dt = (
                + 10**R_MgII_to_MgIII_via_HII(Tx) * nMgII * nHII
                + 10**R_MgII_to_MgIII_via_e(Tx) * nMgII * ne
                - 10**R_MgIII_to_MgII_via_HI(Tx) * nMgIII * nHI
                - 10**R_MgIII_to_MgII_via_e(Tx) * nMgIII * ne
                - 10**R_MgIII_to_MgIV_via_e(Tx) * nMgIII * ne
                + 10**R_MgIV_to_MgIII_via_HeI(Tx) * nMgIV * nHeI
                + 10**R_MgIV_to_MgIII_via_e(Tx) * nMgIV * ne
                + 10**R_MgIV_to_MgIII_via_HI(Tx) * nMgIV * nHI
               )

  dnMgIV_dt = (
               + 10**R_MgIII_to_MgIV_via_e(Tx) * nMgIII * ne
               - 10**R_MgIV_to_MgIII_via_HeI(Tx) * nMgIV * nHeI
               - 10**R_MgIV_to_MgIII_via_e(Tx) * nMgIV * ne
               - 10**R_MgIV_to_MgV_via_e(Tx) * nMgIV * ne
               - 10**R_MgIV_to_MgIII_via_HI(Tx) * nMgIV * nHI
               + 10**R_MgV_to_MgIV_via_HeI(Tx) * nMgV * nHeI
               + 10**R_MgV_to_MgIV_via_HI(Tx) * nMgV * nHI
               + 10**R_MgV_to_MgIV_via_e(Tx) * nMgV * ne
              )

  dnMgV_dt = (
              + 10**R_MgIV_to_MgV_via_e(Tx) * nMgIV * ne
              - 10**R_MgV_to_MgIV_via_HeI(Tx) * nMgV * nHeI
              - 10**R_MgV_to_MgIV_via_HI(Tx) * nMgV * nHI
              - 10**R_MgV_to_MgIV_via_e(Tx) * nMgV * ne
              - 10**R_MgV_to_MgVI_via_e(Tx) * nMgV * ne
              + 10**R_MgVI_to_MgV_via_HI(Tx) * nMgVI * nHI
              + 10**R_MgVI_to_MgV_via_e(Tx) * nMgVI * ne
             )

  dnMgVI_dt = (
               + 10**R_MgV_to_MgVI_via_e(Tx) * nMgV * ne
               - 10**R_MgVI_to_MgV_via_HI(Tx) * nMgVI * nHI
               - 10**R_MgVI_to_MgVII_via_e(Tx) * nMgVI * ne
               - 10**R_MgVI_to_MgV_via_e(Tx) * nMgVI * ne
               + 10**R_MgVII_to_MgVI_via_e(Tx) * nMgVII * ne
              )

  dnMgVII_dt = (
                + 10**R_MgVI_to_MgVII_via_e(Tx) * nMgVI * ne
                - 10**R_MgVII_to_MgVIII_via_e(Tx) * nMgVII * ne
                - 10**R_MgVII_to_MgVI_via_e(Tx) * nMgVII * ne
                + 10**R_MgVIII_to_MgVII_via_e(Tx) * nMgVIII * ne
               )

  dnMgVIII_dt = (
                 + 10**R_MgVII_to_MgVIII_via_e(Tx) * nMgVII * ne
                 - 10**R_MgVIII_to_MgVII_via_e(Tx) * nMgVIII * ne
                 - 10**R_MgVIII_to_MgIX_via_e(Tx) * nMgVIII * ne
                 + 10**R_MgIX_to_MgVIII_via_e(Tx) * nMgIX * ne
                )

  dnMgIX_dt = (
               + 10**R_MgVIII_to_MgIX_via_e(Tx) * nMgVIII * ne
               - 10**R_MgIX_to_MgVIII_via_e(Tx) * nMgIX * ne
               - 10**R_MgIX_to_MgX_via_e(Tx) * nMgIX * ne
               + 10**R_MgX_to_MgIX_via_e(Tx) * nMgX * ne
              )

  dnMgX_dt = (
              + 10**R_MgIX_to_MgX_via_e(Tx) * nMgIX * ne
              - 10**R_MgX_to_MgIX_via_e(Tx) * nMgX * ne
              - 10**R_MgX_to_MgXI_via_e(Tx) * nMgX * ne
              + 10**R_MgXI_to_MgX_via_e(Tx) * nMgXI * ne
             )

  dnMgXI_dt = (
               + 10**R_MgX_to_MgXI_via_e(Tx) * nMgX * ne
               - 10**R_MgXI_to_MgX_via_e(Tx) * nMgXI * ne
               - 10**R_MgXI_to_MgXII_via_e(Tx) * nMgXI * ne
               + 10**R_MgXII_to_MgXI_via_e(Tx) * nMgXII * ne
              )

  dnMgXII_dt = (
                + 10**R_MgXI_to_MgXII_via_e(Tx) * nMgXI * ne
                - 10**R_MgXII_to_MgXIII_via_e(Tx) * nMgXII * ne
                - 10**R_MgXII_to_MgXI_via_e(Tx) * nMgXII * ne
                + 10**R_MgXIII_to_MgXII_via_e(Tx) * nMgXIII * ne
               )

  dnMgXIII_dt = (
                 + 10**R_MgXII_to_MgXIII_via_e(Tx) * nMgXII * ne
                 - 10**R_MgXIII_to_MgXII_via_e(Tx) * nMgXIII * ne
                )

  dnSiI_dt = (
              - 10**R_SiI_to_SiII_via_HeII(Tx) * nSiI * nHeII
              + 10**grain_rec_SiII_to_SiI * nSiII * ne # grain_recombination
              - 10**R_SiI_to_SiII_via_HII(Tx) * nSiI * nHII
              - 10**R_SiI_to_SiII_via_e(Tx) * nSiI * ne
              + 10**R_SiII_to_SiI_via_e(Tx) * nSiII * ne
              - const_CII_SiI_to_CI_SiII * nCII * nSiI # constant rate
              + const_SiII_MgI_to_SiI_MgII * nSiII * nMgI # constant rate
             )

  dnSiII_dt = (
               + 10**R_SiI_to_SiII_via_HeII(Tx) * nSiI * nHeII
               - 10**grain_rec_SiII_to_SiI * nSiII * ne # grain_recombination
               + 10**R_SiI_to_SiII_via_HII(Tx) * nSiI * nHII
               + 10**R_SiI_to_SiII_via_e(Tx) * nSiI * ne
               - 10**R_SiII_to_SiIII_via_HII(Tx) * nSiII * nHII
               - 10**R_SiII_to_SiI_via_e(Tx) * nSiII * ne
               - 10**R_SiII_to_SiIII_via_e(Tx) * nSiII * ne
               - 10**R_SiII_to_SiIII_via_HeII(Tx) * nSiII * nHeII
               + 10**R_SiIII_to_SiII_via_e(Tx) * nSiIII * ne
               + 10**R_SiIII_to_SiII_via_HI(Tx) * nSiIII * nHI
               + const_CII_SiI_to_CI_SiII * nCII * nSiI # constant rate
               - const_SiII_MgI_to_SiI_MgII * nSiII * nMgI # constant rate
              )

  dnSiIII_dt = (
                + 10**R_SiII_to_SiIII_via_HII(Tx) * nSiII * nHII
                + 10**R_SiII_to_SiIII_via_e(Tx) * nSiII * ne
                + 10**R_SiII_to_SiIII_via_HeII(Tx) * nSiII * nHeII
                - 10**R_SiIII_to_SiIV_via_e(Tx) * nSiIII * ne
                - 10**R_SiIII_to_SiII_via_e(Tx) * nSiIII * ne
                - 10**R_SiIII_to_SiII_via_HI(Tx) * nSiIII * nHI
                - 10**R_SiIII_to_SiIV_via_HeII(Tx) * nSiIII * nHeII
                + 10**R_SiIV_to_SiIII_via_HeI(Tx) * nSiIV * nHeI
                + 10**R_SiIV_to_SiIII_via_HI(Tx) * nSiIV * nHI
                + 10**R_SiIV_to_SiIII_via_e(Tx) * nSiIV * ne
               )

  dnSiIV_dt = (
               + 10**R_SiIII_to_SiIV_via_e(Tx) * nSiIII * ne
               + 10**R_SiIII_to_SiIV_via_HeII(Tx) * nSiIII * nHeII
               - 10**R_SiIV_to_SiIII_via_HeI(Tx) * nSiIV * nHeI
               - 10**R_SiIV_to_SiIII_via_HI(Tx) * nSiIV * nHI
               - 10**R_SiIV_to_SiIII_via_e(Tx) * nSiIV * ne
               - 10**R_SiIV_to_SiV_via_e(Tx) * nSiIV * ne
               + 10**R_SiV_to_SiIV_via_HI(Tx) * nSiV * nHI
               + 10**R_SiV_to_SiIV_via_e(Tx) * nSiV * ne
               + 10**R_SiV_to_SiIV_via_HeI(Tx) * nSiV * nHeI
              )

  dnSiV_dt = (
              + 10**R_SiIV_to_SiV_via_e(Tx) * nSiIV * ne
              - 10**R_SiV_to_SiIV_via_HI(Tx) * nSiV * nHI
              - 10**R_SiV_to_SiIV_via_e(Tx) * nSiV * ne
              - 10**R_SiV_to_SiVI_via_e(Tx) * nSiV * ne
              - 10**R_SiV_to_SiIV_via_HeI(Tx) * nSiV * nHeI
              + 10**R_SiVI_to_SiV_via_HI(Tx) * nSiVI * nHI
              + 10**R_SiVI_to_SiV_via_e(Tx) * nSiVI * ne
             )

  dnSiVI_dt = (
               + 10**R_SiV_to_SiVI_via_e(Tx) * nSiV * ne
               - 10**R_SiVI_to_SiV_via_HI(Tx) * nSiVI * nHI
               - 10**R_SiVI_to_SiV_via_e(Tx) * nSiVI * ne
               - 10**R_SiVI_to_SiVII_via_e(Tx) * nSiVI * ne
               + 10**R_SiVII_to_SiVI_via_e(Tx) * nSiVII * ne
              )

  dnSiVII_dt = (
                + 10**R_SiVI_to_SiVII_via_e(Tx) * nSiVI * ne
                - 10**R_SiVII_to_SiVI_via_e(Tx) * nSiVII * ne
                - 10**R_SiVII_to_SiVIII_via_e(Tx) * nSiVII * ne
                + 10**R_SiVIII_to_SiVII_via_e(Tx) * nSiVIII * ne
               )

  dnSiVIII_dt = (
                 + 10**R_SiVII_to_SiVIII_via_e(Tx) * nSiVII * ne
                 - 10**R_SiVIII_to_SiVII_via_e(Tx) * nSiVIII * ne
                 - 10**R_SiVIII_to_SiIX_via_e(Tx) * nSiVIII * ne
                 + 10**R_SiIX_to_SiVIII_via_e(Tx) * nSiIX * ne
                )

  dnSiIX_dt = (
               + 10**R_SiVIII_to_SiIX_via_e(Tx) * nSiVIII * ne
               - 10**R_SiIX_to_SiX_via_e(Tx) * nSiIX * ne
               - 10**R_SiIX_to_SiVIII_via_e(Tx) * nSiIX * ne
               + 10**R_SiX_to_SiIX_via_e(Tx) * nSiX * ne
              )

  dnSiX_dt = (
              + 10**R_SiIX_to_SiX_via_e(Tx) * nSiIX * ne
              - 10**R_SiX_to_SiXI_via_e(Tx) * nSiX * ne
              - 10**R_SiX_to_SiIX_via_e(Tx) * nSiX * ne
              + 10**R_SiXI_to_SiX_via_e(Tx) * nSiXI * ne
             )

  dnSiXI_dt = (
               + 10**R_SiX_to_SiXI_via_e(Tx) * nSiX * ne
               - 10**R_SiXI_to_SiX_via_e(Tx) * nSiXI * ne
               - 10**R_SiXI_to_SiXII_via_e(Tx) * nSiXI * ne
               + 10**R_SiXII_to_SiXI_via_e(Tx) * nSiXII * ne
              )

  dnSiXII_dt = (
                + 10**R_SiXI_to_SiXII_via_e(Tx) * nSiXI * ne
                - 10**R_SiXII_to_SiXI_via_e(Tx) * nSiXII * ne
                - 10**R_SiXII_to_SiXIII_via_e(Tx) * nSiXII * ne
                + 10**R_SiXIII_to_SiXII_via_e(Tx) * nSiXIII * ne
               )

  dnSiXIII_dt = (
                 + 10**R_SiXII_to_SiXIII_via_e(Tx) * nSiXII * ne
                 - 10**R_SiXIII_to_SiXII_via_e(Tx) * nSiXIII * ne
                 - 10**R_SiXIII_to_SiXIV_via_e(Tx) * nSiXIII * ne
                 + 10**R_SiXIV_to_SiXIII_via_e(Tx) * nSiXIV * ne
                )

  dnSiXIV_dt = (
                + 10**R_SiXIII_to_SiXIV_via_e(Tx) * nSiXIII * ne
                - 10**R_SiXIV_to_SiXIII_via_e(Tx) * nSiXIV * ne
                - 10**R_SiXIV_to_SiXV_via_e(Tx) * nSiXIV * ne
                + 10**R_SiXV_to_SiXIV_via_e(Tx) * nSiXV * ne
               )

  dnSiXV_dt = (
               + 10**R_SiXIV_to_SiXV_via_e(Tx) * nSiXIV * ne
               - 10**R_SiXV_to_SiXIV_via_e(Tx) * nSiXV * ne
              )

  dnSI_dt = (
             - 10**R_SI_to_SII_via_HII(Tx) * nSI * nHII
             + 10**grain_rec_SII_to_SI * nSII * ne # grain_recombination
             - 10**R_SI_to_SII_via_e(Tx) * nSI * ne
             + 10**R_SII_to_SI_via_HI(Tx) * nSII * nHI
             + 10**R_SII_to_SI_via_e(Tx) * nSII * ne
             + const_SII_MgI_to_SI_MgII * nSII * nMgI # constant rate
            )

  dnSII_dt = (
              + 10**R_SI_to_SII_via_HII(Tx) * nSI * nHII
              - 10**grain_rec_SII_to_SI * nSII * ne # grain_recombination
              + 10**R_SI_to_SII_via_e(Tx) * nSI * ne
              - 10**R_SII_to_SI_via_HI(Tx) * nSII * nHI
              - 10**R_SII_to_SIII_via_HeII(Tx) * nSII * nHeII
              - 10**R_SII_to_SI_via_e(Tx) * nSII * ne
              - 10**R_SII_to_SIII_via_e(Tx) * nSII * ne
              + 10**R_SIII_to_SII_via_e(Tx) * nSIII * ne
              + 10**R_SIII_to_SII_via_HI(Tx) * nSIII * nHI
              - const_SII_MgI_to_SI_MgII * nSII * nMgI # constant rate
             )

  dnSIII_dt = (
               + 10**R_SII_to_SIII_via_HeII(Tx) * nSII * nHeII
               + 10**R_SII_to_SIII_via_e(Tx) * nSII * ne
               - 10**R_SIII_to_SIV_via_HeII(Tx) * nSIII * nHeII
               - 10**R_SIII_to_SII_via_e(Tx) * nSIII * ne
               - 10**R_SIII_to_SIV_via_e(Tx) * nSIII * ne
               - 10**R_SIII_to_SII_via_HI(Tx) * nSIII * nHI
               + 10**R_SIV_to_SIII_via_e(Tx) * nSIV * ne
               + 10**R_SIV_to_SIII_via_HI(Tx) * nSIV * nHI
               + 10**R_SIV_to_SIII_via_HeI(Tx) * nSIV * nHeI
              )

  dnSIV_dt = (
              + 10**R_SIII_to_SIV_via_HeII(Tx) * nSIII * nHeII
              + 10**R_SIII_to_SIV_via_e(Tx) * nSIII * ne
              - 10**R_SIV_to_SV_via_e(Tx) * nSIV * ne
              - 10**R_SIV_to_SIII_via_e(Tx) * nSIV * ne
              - 10**R_SIV_to_SIII_via_HI(Tx) * nSIV * nHI
              - 10**R_SIV_to_SIII_via_HeI(Tx) * nSIV * nHeI
              + 10**R_SV_to_SIV_via_HeI(Tx) * nSV * nHeI
              + 10**R_SV_to_SIV_via_HI(Tx) * nSV * nHI
              + 10**R_SV_to_SIV_via_e(Tx) * nSV * ne
             )

  dnSV_dt = (
             + 10**R_SIV_to_SV_via_e(Tx) * nSIV * ne
             - 10**R_SV_to_SIV_via_HeI(Tx) * nSV * nHeI
             - 10**R_SV_to_SIV_via_HI(Tx) * nSV * nHI
             - 10**R_SV_to_SIV_via_e(Tx) * nSV * ne
             - 10**R_SV_to_SVI_via_e(Tx) * nSV * ne
             + 10**R_SVI_to_SV_via_HI(Tx) * nSVI * nHI
             + 10**R_SVI_to_SV_via_e(Tx) * nSVI * ne
            )

  dnSVI_dt = (
              + 10**R_SV_to_SVI_via_e(Tx) * nSV * ne
              - 10**R_SVI_to_SV_via_HI(Tx) * nSVI * nHI
              - 10**R_SVI_to_SV_via_e(Tx) * nSVI * ne
              - 10**R_SVI_to_SVII_via_e(Tx) * nSVI * ne
              + 10**R_SVII_to_SVI_via_e(Tx) * nSVII * ne
             )

  dnSVII_dt = (
               + 10**R_SVI_to_SVII_via_e(Tx) * nSVI * ne
               - 10**R_SVII_to_SVI_via_e(Tx) * nSVII * ne
               - 10**R_SVII_to_SVIII_via_e(Tx) * nSVII * ne
               + 10**R_SVIII_to_SVII_via_e(Tx) * nSVIII * ne
              )

  dnSVIII_dt = (
                + 10**R_SVII_to_SVIII_via_e(Tx) * nSVII * ne
                - 10**R_SVIII_to_SVII_via_e(Tx) * nSVIII * ne
                - 10**R_SVIII_to_SIX_via_e(Tx) * nSVIII * ne
                + 10**R_SIX_to_SVIII_via_e(Tx) * nSIX * ne
               )

  dnSIX_dt = (
              + 10**R_SVIII_to_SIX_via_e(Tx) * nSVIII * ne
              - 10**R_SIX_to_SVIII_via_e(Tx) * nSIX * ne
              - 10**R_SIX_to_SX_via_e(Tx) * nSIX * ne
              + 10**R_SX_to_SIX_via_e(Tx) * nSX * ne
             )

  dnSX_dt = (
             + 10**R_SIX_to_SX_via_e(Tx) * nSIX * ne
             - 10**R_SX_to_SXI_via_e(Tx) * nSX * ne
             - 10**R_SX_to_SIX_via_e(Tx) * nSX * ne
             + 10**R_SXI_to_SX_via_e(Tx) * nSXI * ne
            )

  dnSXI_dt = (
              + 10**R_SX_to_SXI_via_e(Tx) * nSX * ne
              - 10**R_SXI_to_SX_via_e(Tx) * nSXI * ne
              - 10**R_SXI_to_SXII_via_e(Tx) * nSXI * ne
              + 10**R_SXII_to_SXI_via_e(Tx) * nSXII * ne
             )

  dnSXII_dt = (
               + 10**R_SXI_to_SXII_via_e(Tx) * nSXI * ne
               - 10**R_SXII_to_SXIII_via_e(Tx) * nSXII * ne
               - 10**R_SXII_to_SXI_via_e(Tx) * nSXII * ne
               + 10**R_SXIII_to_SXII_via_e(Tx) * nSXIII * ne
              )

  dnSXIII_dt = (
                + 10**R_SXII_to_SXIII_via_e(Tx) * nSXII * ne
                - 10**R_SXIII_to_SXII_via_e(Tx) * nSXIII * ne
                - 10**R_SXIII_to_SXIV_via_e(Tx) * nSXIII * ne
                + 10**R_SXIV_to_SXIII_via_e(Tx) * nSXIV * ne
               )

  dnSXIV_dt = (
               + 10**R_SXIII_to_SXIV_via_e(Tx) * nSXIII * ne
               - 10**R_SXIV_to_SXIII_via_e(Tx) * nSXIV * ne
               - 10**R_SXIV_to_SXV_via_e(Tx) * nSXIV * ne
               + 10**R_SXV_to_SXIV_via_e(Tx) * nSXV * ne
              )

  dnSXV_dt = (
              + 10**R_SXIV_to_SXV_via_e(Tx) * nSXIV * ne
              - 10**R_SXV_to_SXIV_via_e(Tx) * nSXV * ne
              - 10**R_SXV_to_SXVI_via_e(Tx) * nSXV * ne
              + 10**R_SXVI_to_SXV_via_e(Tx) * nSXVI * ne
             )

  dnSXVI_dt = (
               + 10**R_SXV_to_SXVI_via_e(Tx) * nSXV * ne
               - 10**R_SXVI_to_SXV_via_e(Tx) * nSXVI * ne
               - 10**R_SXVI_to_SXVII_via_e(Tx) * nSXVI * ne
               + 10**R_SXVII_to_SXVI_via_e(Tx) * nSXVII * ne
              )

  dnSXVII_dt = (
                + 10**R_SXVI_to_SXVII_via_e(Tx) * nSXVI * ne
                - 10**R_SXVII_to_SXVI_via_e(Tx) * nSXVII * ne
               )

  Lamb = Lambda(
                T, nHI, nHII, nHm, nHeI, nHeII, nHeIII, nCI, nCII, nCIII, nCIV, 
                nCV, nCVI, nCVII, nCm, nNI, nNII, nNIII, nNIV, nNV, nNVI, 
                nNVII, nNVIII, nOI, nOII, nOIII, nOIV, nOV, nOVI, nOVII, nOVIII, 
                nOIX, nOm, nNeI, nNeII, nNeIII, nNeIV, nNeV, nNeVI, nNeVII, nNeVIII, 
                nNeIX, nNeX, nNeXI, nMgI, nMgII, nMgIII, nMgIV, nMgV, nMgVI, nMgVII, 
                nMgVIII, nMgIX, nMgX, nMgXI, nMgXII, nMgXIII, nSiI, nSiII, nSiIII, nSiIV, 
                nSiV, nSiVI, nSiVII, nSiVIII, nSiIX, nSiX, nSiXI, nSiXII, nSiXIII, nSiXIV, 
                nSiXV, nSI, nSII, nSIII, nSIV, nSV, nSVI, nSVII, nSVIII, nSIX, 
                nSX, nSXI, nSXII, nSXIII, nSXIV, nSXV, nSXVI, nSXVII)

  dne_dt = (
           + 1 * dnHII_dt + -1 * dnHm_dt + 1 * dnHeII_dt + 2 * dnHeIII_dt + 1 * dnCII_dt + 2 * dnCIII_dt
           + 3 * dnCIV_dt + 4 * dnCV_dt + 5 * dnCVI_dt + 6 * dnCVII_dt + -1 * dnCm_dt + 1 * dnNII_dt
           + 2 * dnNIII_dt + 3 * dnNIV_dt + 4 * dnNV_dt + 5 * dnNVI_dt + 6 * dnNVII_dt + 7 * dnNVIII_dt
           + 1 * dnOII_dt + 2 * dnOIII_dt + 3 * dnOIV_dt + 4 * dnOV_dt + 5 * dnOVI_dt + 6 * dnOVII_dt
           + 7 * dnOVIII_dt + 8 * dnOIX_dt + -1 * dnOm_dt + 1 * dnNeII_dt + 2 * dnNeIII_dt + 3 * dnNeIV_dt
           + 4 * dnNeV_dt + 5 * dnNeVI_dt + 6 * dnNeVII_dt + 7 * dnNeVIII_dt + 8 * dnNeIX_dt + 9 * dnNeX_dt
           + 10 * dnNeXI_dt + 1 * dnMgII_dt + 2 * dnMgIII_dt + 3 * dnMgIV_dt + 4 * dnMgV_dt + 5 * dnMgVI_dt
           + 6 * dnMgVII_dt + 7 * dnMgVIII_dt + 8 * dnMgIX_dt + 9 * dnMgX_dt + 10 * dnMgXI_dt + 11 * dnMgXII_dt
           + 12 * dnMgXIII_dt + 1 * dnSiII_dt + 2 * dnSiIII_dt + 3 * dnSiIV_dt + 4 * dnSiV_dt + 5 * dnSiVI_dt
           + 6 * dnSiVII_dt + 7 * dnSiVIII_dt + 8 * dnSiIX_dt + 9 * dnSiX_dt + 10 * dnSiXI_dt + 11 * dnSiXII_dt
           + 12 * dnSiXIII_dt + 13 * dnSiXIV_dt + 14 * dnSiXV_dt + 1 * dnSII_dt + 2 * dnSIII_dt + 3 * dnSIV_dt
           + 4 * dnSV_dt + 5 * dnSVI_dt + 6 * dnSVII_dt + 7 * dnSVIII_dt + 8 * dnSIX_dt + 9 * dnSX_dt
           + 10 * dnSXI_dt + 11 * dnSXII_dt + 12 * dnSXIII_dt + 13 * dnSXIV_dt + 14 * dnSXV_dt + 15 * dnSXVI_dt
           + 16 * dnSXVII_dt
           )

  dntot_dt = (
                dne_dt + dnHI_dt + dnHII_dt + dnHm_dt + dnHeI_dt + dnHeII_dt + dnHeIII_dt
              + dnCI_dt + dnCII_dt + dnCIII_dt + dnCIV_dt + dnCV_dt + dnCVI_dt
              + dnCVII_dt + dnCm_dt + dnNI_dt + dnNII_dt + dnNIII_dt + dnNIV_dt
              + dnNV_dt + dnNVI_dt + dnNVII_dt + dnNVIII_dt + dnOI_dt + dnOII_dt
              + dnOIII_dt + dnOIV_dt + dnOV_dt + dnOVI_dt + dnOVII_dt + dnOVIII_dt
              + dnOIX_dt + dnOm_dt + dnNeI_dt + dnNeII_dt + dnNeIII_dt + dnNeIV_dt
              + dnNeV_dt + dnNeVI_dt + dnNeVII_dt + dnNeVIII_dt + dnNeIX_dt + dnNeX_dt
              + dnNeXI_dt + dnMgI_dt + dnMgII_dt + dnMgIII_dt + dnMgIV_dt + dnMgV_dt
              + dnMgVI_dt + dnMgVII_dt + dnMgVIII_dt + dnMgIX_dt + dnMgX_dt + dnMgXI_dt
              + dnMgXII_dt + dnMgXIII_dt + dnSiI_dt + dnSiII_dt + dnSiIII_dt + dnSiIV_dt
              + dnSiV_dt + dnSiVI_dt + dnSiVII_dt + dnSiVIII_dt + dnSiIX_dt + dnSiX_dt
              + dnSiXI_dt + dnSiXII_dt + dnSiXIII_dt + dnSiXIV_dt + dnSiXV_dt + dnSI_dt
              + dnSII_dt + dnSIII_dt + dnSIV_dt + dnSV_dt + dnSVI_dt + dnSVII_dt
              + dnSVIII_dt + dnSIX_dt + dnSX_dt + dnSXI_dt + dnSXII_dt + dnSXIII_dt
              + dnSXIV_dt + dnSXV_dt + dnSXVI_dt + dnSXVII_dt
             )

  dT_dt = -1.0 * (gamma - 1.0) / kB / ntot * (Lamb + 1. / (gamma - 1.) * kB * T * dntot_dt)

  return [
          dnHI_dt, dnHII_dt, dnHm_dt, dnHeI_dt, dnHeII_dt, dnHeIII_dt,
          dnCI_dt, dnCII_dt, dnCIII_dt, dnCIV_dt, dnCV_dt, dnCVI_dt,
          dnCVII_dt, dnCm_dt, dnNI_dt, dnNII_dt, dnNIII_dt, dnNIV_dt,
          dnNV_dt, dnNVI_dt, dnNVII_dt, dnNVIII_dt, dnOI_dt, dnOII_dt,
          dnOIII_dt, dnOIV_dt, dnOV_dt, dnOVI_dt, dnOVII_dt, dnOVIII_dt,
          dnOIX_dt, dnOm_dt, dnNeI_dt, dnNeII_dt, dnNeIII_dt, dnNeIV_dt,
          dnNeV_dt, dnNeVI_dt, dnNeVII_dt, dnNeVIII_dt, dnNeIX_dt, dnNeX_dt,
          dnNeXI_dt, dnMgI_dt, dnMgII_dt, dnMgIII_dt, dnMgIV_dt, dnMgV_dt,
          dnMgVI_dt, dnMgVII_dt, dnMgVIII_dt, dnMgIX_dt, dnMgX_dt, dnMgXI_dt,
          dnMgXII_dt, dnMgXIII_dt, dnSiI_dt, dnSiII_dt, dnSiIII_dt, dnSiIV_dt,
          dnSiV_dt, dnSiVI_dt, dnSiVII_dt, dnSiVIII_dt, dnSiIX_dt, dnSiX_dt,
          dnSiXI_dt, dnSiXII_dt, dnSiXIII_dt, dnSiXIV_dt, dnSiXV_dt, dnSI_dt,
          dnSII_dt, dnSIII_dt, dnSIV_dt, dnSV_dt, dnSVI_dt, dnSVII_dt,
          dnSVIII_dt, dnSIX_dt, dnSX_dt, dnSXI_dt, dnSXII_dt, dnSXIII_dt,
          dnSXIV_dt, dnSXV_dt, dnSXVI_dt, dnSXVII_dt, dT_dt
         ]








TA = time.time()

print("Running ...")

nH = 1000.0

He_solar = 10**(-1.00)
nHe = He_solar * nH

C_solar = 10**(-3.61)
nC = C_solar * nH

N_solar = 10**(-4.07)
nN = N_solar * nH

O_solar = 10**(-3.31)
nO = O_solar * nH

Ne_solar = 10**(-4.00)
nNe = Ne_solar * nH

Mg_solar = 10**(-4.46)
nMg = Mg_solar * nH

Si_solar = 10**(-4.46)
nSi = Si_solar * nH

S_solar = 10**(-4.73)
nS = S_solar * nH

T_i = 10**7.00

nHm_i = 1e-5 * nH
nH0_i = 0.001 * nH
nHp_i = nH - nH0_i

nHe0_i = 0.0001 * nHe
nHep_i = 0.001 * nHe
nHepp_i= nHe - nHe0_i - nHep_i

nC0_i = 1e-6 * nC
nC1_i = 1e-5 * nC
nC2_i = 1e-4 * nC
nC3_i = 1e-3 * nC
nC4_i = 1e-2 * nC
nC5_i = 1e-2 * nC
nCm_i = 1e-6 * nC
nC6_i = nC - (nC0_i + nC1_i + nC2_i + nC3_i + nC4_i + nC5_i + nCm_i)

nN0_i = 1e-7 * nN
nN1_i = 1e-6 * nN
nN2_i = 1e-5 * nN
nN3_i = 1e-4 * nN
nN4_i = 1e-3 * nN
nN5_i = 1e-2 * nN
nN6_i = 1e-2 * nN
nN7_i = nN - (nN0_i + nN1_i + nN2_i + nN3_i + nN4_i + nN5_i + nN6_i)

nO0_i = 1e-8 * nO
nO1_i = 1e-7 * nO
nO2_i = 1e-6 * nO
nO3_i = 1e-5 * nO
nO4_i = 1e-4 * nO
nO5_i = 1e-3 * nO
nO6_i = 1e-2 * nO
nO7_i = 1e-2 * nO
nOm_i = 1e-6 * nO
nO8_i = nO - (nO0_i + nO1_i + nO2_i + nO3_i + nO4_i + nO5_i + nO6_i + nO7_i + nOm_i)

nNe0_i = 1e-10 * nNe
nNe1_i = 1e-9 * nNe
nNe2_i = 1e-8 * nNe
nNe3_i = 1e-7 * nNe
nNe4_i = 1e-6 * nNe
nNe5_i = 1e-5 * nNe
nNe6_i = 1e-4 * nNe
nNe7_i = 1e-3 * nNe
nNe8_i = 1e-2 * nNe
nNe9_i = 1e-2 * nNe
nNe10_i = nNe - (nNe0_i + nNe1_i + nNe2_i + nNe3_i + nNe4_i + nNe5_i + nNe6_i + nNe7_i + nNe8_i + nNe9_i)

nMg0_i = 1e-12 * nMg
nMg1_i = 1e-11 * nMg
nMg2_i = 1e-10 * nMg
nMg3_i = 1e-9 * nMg
nMg4_i = 1e-8 * nMg
nMg5_i = 1e-7 * nMg
nMg6_i = 1e-6 * nMg
nMg7_i = 1e-5 * nMg
nMg8_i = 1e-4 * nMg
nMg9_i = 1e-3 * nMg
nMg10_i = 1e-2 * nMg
nMg11_i = 1e-2 * nMg
nMg12_i = nMg - (nMg0_i + nMg1_i + nMg2_i + nMg3_i + nMg4_i + nMg5_i + nMg6_i + nMg7_i + nMg8_i + nMg9_i + nMg10_i + nMg11_i)

nSi0_i = 1e-14 * nSi
nSi1_i = 1e-13 * nSi
nSi2_i = 1e-12 * nSi
nSi3_i = 1e-11 * nSi
nSi4_i = 1e-10 * nSi
nSi5_i = 1e-9 * nSi
nSi6_i = 1e-8 * nSi
nSi7_i = 1e-7 * nSi
nSi8_i = 1e-6 * nSi
nSi9_i = 1e-5 * nSi
nSi10_i = 1e-4 * nSi
nSi11_i = 1e-3 * nSi
nSi12_i = 1e-2 * nSi
nSi13_i = 1e-2 * nSi
nSi14_i = nSi - (nSi0_i + nSi1_i + nSi2_i + nSi3_i + nSi4_i + nSi5_i + nSi6_i + nSi7_i + nSi8_i + nSi9_i + nSi10_i + nSi11_i + nSi12_i + nSi13_i)

nS0_i = 1e-16 * nS
nS1_i = 1e-15 * nS
nS2_i = 1e-14 * nS
nS3_i = 1e-13 * nS
nS4_i = 1e-12 * nS
nS5_i = 1e-11 * nS
nS6_i = 1e-10 * nS
nS7_i = 1e-9 * nS
nS8_i = 1e-8 * nS
nS9_i = 1e-7 * nS
nS10_i = 1e-6 * nS
nS11_i = 1e-5 * nS
nS12_i = 1e-4 * nS
nS13_i = 1e-3 * nS
nS14_i = 1e-2 * nS
nS15_i = 1e-2 * nS
nS16_i = nS - (nS0_i + nS1_i + nS2_i + nS3_i + nS4_i + nS5_i + nS6_i + nS7_i + nS8_i + nS9_i + nS10_i + nS11_i + nS12_i + nS13_i + nS14_i + nS15_i)

y0 = [
      nH0_i, nHp_i, nHm_i, nHe0_i, nHep_i, nHepp_i, 
      nC0_i, nC1_i, nC2_i, nC3_i, nC4_i, nC5_i, 
      nC6_i, nCm_i, nN0_i, nN1_i, nN2_i, nN3_i, 
      nN4_i, nN5_i, nN6_i, nN7_i, nO0_i, nO1_i, 
      nO2_i, nO3_i, nO4_i, nO5_i, nO6_i, nO7_i, 
      nO8_i, nOm_i, nNe0_i, nNe1_i, nNe2_i, nNe3_i, 
      nNe4_i, nNe5_i, nNe6_i, nNe7_i, nNe8_i, nNe9_i, 
      nNe10_i, nMg0_i, nMg1_i, nMg2_i, nMg3_i, nMg4_i, 
      nMg5_i, nMg6_i, nMg7_i, nMg8_i, nMg9_i, nMg10_i, 
      nMg11_i, nMg12_i, nSi0_i, nSi1_i, nSi2_i, nSi3_i, 
      nSi4_i, nSi5_i, nSi6_i, nSi7_i, nSi8_i, nSi9_i, 
      nSi10_i, nSi11_i, nSi12_i, nSi13_i, nSi14_i, nS0_i, 
      nS1_i, nS2_i, nS3_i, nS4_i, nS5_i, nS6_i, 
      nS7_i, nS8_i, nS9_i, nS10_i, nS11_i, nS12_i, 
      nS13_i, nS14_i, nS15_i, nS16_i, 
      T_i
      ]

A_v = 1.0
G0 = 0.01
dust_ratio = 0.01

t_span = (1*3.16e7, 8000*3.16e7)

solution = solve_ivp(func, t_span, y0, method="LSODA", dense_output=True)

t = np.linspace(t_span[0], t_span[1], 50000) # This 10000 is not years, it is the number of points in linspace !!!!
y = solution.sol(t)

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


nO0 = y[22, :]
nO1 = y[23, :]
nO2 = y[24, :]
nO3 = y[25, :]
nO4 = y[26, :]
nO5 = y[27, :]
nO6 = y[28, :]
nO7 = y[29, :]
nO8 = y[30, :]
nOm = y[31, :]


nNe0 = y[32, :]
nNe1 = y[33, :]
nNe2 = y[34, :]
nNe3 = y[35, :]
nNe4 = y[36, :]
nNe5 = y[37, :]
nNe6 = y[38, :]
nNe7 = y[39, :]
nNe8 = y[40, :]
nNe9 = y[41, :]
nNe10 = y[42, :]


nMg0 = y[43, :]
nMg1 = y[44, :]
nMg2 = y[45, :]
nMg3 = y[46, :]
nMg4 = y[47, :]
nMg5 = y[48, :]
nMg6 = y[49, :]
nMg7 = y[50, :]
nMg8 = y[51, :]
nMg9 = y[52, :]
nMg10 = y[53, :]
nMg11 = y[54, :]
nMg12 = y[55, :]


nSi0 = y[56, :]
nSi1 = y[57, :]
nSi2 = y[58, :]
nSi3 = y[59, :]
nSi4 = y[60, :]
nSi5 = y[61, :]
nSi6 = y[62, :]
nSi7 = y[63, :]
nSi8 = y[64, :]
nSi9 = y[65, :]
nSi10 = y[66, :]
nSi11 = y[67, :]
nSi12 = y[68, :]
nSi13 = y[69, :]
nSi14 = y[70, :]


nS0 = y[71, :]
nS1 = y[72, :]
nS2 = y[73, :]
nS3 = y[74, :]
nS4 = y[75, :]
nS5 = y[76, :]
nS6 = y[77, :]
nS7 = y[78, :]
nS8 = y[79, :]
nS9 = y[80, :]
nS10 = y[81, :]
nS11 = y[82, :]
nS12 = y[83, :]
nS13 = y[84, :]
nS14 = y[85, :]
nS15 = y[86, :]
nS16 = y[87, :]


T = y[88, :]

with open ("chimesData.pkl", "rb") as f:
  data = pickle.load(f)

TEvolx = data["TempEvol"]
AbundEvol = data["chimesAbundEvol"]
t_Arr_in_yrsx = data["t_Arr_in_yrs"]

nH0x   = AbundEvol[1, :]
nHpx   = AbundEvol[2, :]
nHe0x  = AbundEvol[4, :]
nHepx  = AbundEvol[5, :]
nHeppx = AbundEvol[6, :]

nC0x = AbundEvol[7, :]
nC1x = AbundEvol[8, :]
nC2x = AbundEvol[9, :]
nC3x = AbundEvol[10, :]
nC4x = AbundEvol[11, :]
nC5x = AbundEvol[12, :]
nC6x = AbundEvol[13, :]
nCmx = AbundEvol[14, :]
nCx = nC0x + nC1x + nC2x + nC3x + nC4x + nC5x + nC6x + nCmx 

nN0x = AbundEvol[15, :]
nN1x = AbundEvol[16, :]
nN2x = AbundEvol[17, :]
nN3x = AbundEvol[18, :]
nN4x = AbundEvol[19, :]
nN5x = AbundEvol[20, :]
nN6x = AbundEvol[21, :]
nN7x = AbundEvol[22, :]
nNx = nN0x + nN1x + nN2x + nN3x + nN4x + nN5x + nN6x + nN7x 

nO0x = AbundEvol[23, :]
nO1x = AbundEvol[24, :]
nO2x = AbundEvol[25, :]
nO3x = AbundEvol[26, :]
nO4x = AbundEvol[27, :]
nO5x = AbundEvol[28, :]
nO6x = AbundEvol[29, :]
nO7x = AbundEvol[30, :]
nO8x = AbundEvol[31, :]
nOmx = AbundEvol[32, :]
nOx = nO0x + nO1x + nO2x + nO3x + nO4x + nO5x + nO6x + nO7x + nO8x + nOmx 

nNe0x = AbundEvol[33, :]
nNe1x = AbundEvol[34, :]
nNe2x = AbundEvol[35, :]
nNe3x = AbundEvol[36, :]
nNe4x = AbundEvol[37, :]
nNe5x = AbundEvol[38, :]
nNe6x = AbundEvol[39, :]
nNe7x = AbundEvol[40, :]
nNe8x = AbundEvol[41, :]
nNe9x = AbundEvol[42, :]
nNe10x = AbundEvol[43, :]
nNex = nNe0x + nNe1x + nNe2x + nNe3x + nNe4x + nNe5x + nNe6x + nNe7x + nNe8x + nNe9x + nNe10x 

nMg0x = AbundEvol[44, :]
nMg1x = AbundEvol[45, :]
nMg2x = AbundEvol[46, :]
nMg3x = AbundEvol[47, :]
nMg4x = AbundEvol[48, :]
nMg5x = AbundEvol[49, :]
nMg6x = AbundEvol[50, :]
nMg7x = AbundEvol[51, :]
nMg8x = AbundEvol[52, :]
nMg9x = AbundEvol[53, :]
nMg10x = AbundEvol[54, :]
nMg11x = AbundEvol[55, :]
nMg12x = AbundEvol[56, :]
nMgx = nMg0x + nMg1x + nMg2x + nMg3x + nMg4x + nMg5x + nMg6x + nMg7x + nMg8x + nMg9x + nMg10x + nMg11x + nMg12x 

nSi0x = AbundEvol[57, :]
nSi1x = AbundEvol[58, :]
nSi2x = AbundEvol[59, :]
nSi3x = AbundEvol[60, :]
nSi4x = AbundEvol[61, :]
nSi5x = AbundEvol[62, :]
nSi6x = AbundEvol[63, :]
nSi7x = AbundEvol[64, :]
nSi8x = AbundEvol[65, :]
nSi9x = AbundEvol[66, :]
nSi10x = AbundEvol[67, :]
nSi11x = AbundEvol[68, :]
nSi12x = AbundEvol[69, :]
nSi13x = AbundEvol[70, :]
nSi14x = AbundEvol[71, :]
nSix = nSi0x + nSi1x + nSi2x + nSi3x + nSi4x + nSi5x + nSi6x + nSi7x + nSi8x + nSi9x + nSi10x + nSi11x + nSi12x + nSi13x + nSi14x 

nS0x = AbundEvol[72, :]
nS1x = AbundEvol[73, :]
nS2x = AbundEvol[74, :]
nS3x = AbundEvol[75, :]
nS4x = AbundEvol[76, :]
nS5x = AbundEvol[77, :]
nS6x = AbundEvol[78, :]
nS7x = AbundEvol[79, :]
nS8x = AbundEvol[80, :]
nS9x = AbundEvol[81, :]
nS10x = AbundEvol[82, :]
nS11x = AbundEvol[83, :]
nS12x = AbundEvol[84, :]
nS13x = AbundEvol[85, :]
nS14x = AbundEvol[86, :]
nS15x = AbundEvol[87, :]
nS16x = AbundEvol[88, :]
nSx = nS0x + nS1x + nS2x + nS3x + nS4x + nS5x + nS6x + nS7x + nS8x + nS9x + nS10x + nS11x + nS12x + nS13x + nS14x + nS15x + nS16x 

plt.figure(figsize=(10, 5))
plt.scatter(t_yrs, np.log10(T), s=2, color="k", label="my code")
plt.scatter(t_Arr_in_yrsx, np.log10(TEvolx), s=2, color="orange", label="chimes result", linestyle="--")
plt.xlim(0, 10000)
plt.ylim(1, 8)
plt.legend()
plt.savefig("Temp_vs_time.png", dpi = 300)
plt.close()

plt.plot(t_yrs, nHe0, color = 'r', label = 'nHe0')
plt.plot(t_yrs, nHep, color = 'g', label = 'nHep')
plt.plot(t_yrs, nHepp, color = 'b', label = 'nHepp')
plt.plot(t_Arr_in_yrsx, nHe0x, color = 'r', label = 'nHe0 - chimes', linestyle = ':')
plt.plot(t_Arr_in_yrsx, nHepx, color = 'g', label = 'nHep - chimes', linestyle = ':')
plt.plot(t_Arr_in_yrsx, nHeppx, color = 'b', label = 'nHepp - chimes', linestyle = ':')
plt.xlim(0, 10000)
plt.ylim(1e-8, 300)
plt.yscale('log')
plt.title('solve_ivp')
plt.legend()
plt.savefig("nH_He_vs_time.png", dpi = 300)
plt.close()

plt.plot(T, nC0/nC, label = 'nC0', color = '#1f77b4')
plt.plot(TEvolx, nC0x/nCx, label = 'nC0x', color = '#1f77b4', linestyle = ':')
plt.plot(T, nC1/nC, label = 'nC1', color = '#ff7f0e')
plt.plot(TEvolx, nC1x/nCx, label = 'nC1x', color = '#ff7f0e', linestyle = ':')
plt.plot(T, nC2/nC, label = 'nC2', color = '#2ca02c')
plt.plot(TEvolx, nC2x/nCx, label = 'nC2x', color = '#2ca02c', linestyle = ':')
plt.plot(T, nC3/nC, label = 'nC3', color = '#d62728')
plt.plot(TEvolx, nC3x/nCx, label = 'nC3x', color = '#d62728', linestyle = ':')
plt.plot(T, nC4/nC, label = 'nC4', color = '#9467bd')
plt.plot(TEvolx, nC4x/nCx, label = 'nC4x', color = '#9467bd', linestyle = ':')
plt.plot(T, nC5/nC, label = 'nC5', color = '#8c564b')
plt.plot(TEvolx, nC5x/nCx, label = 'nC5x', color = '#8c564b', linestyle = ':')
plt.plot(T, nC6/nC, label = 'nC6', color = '#e377c2')
plt.plot(TEvolx, nC6x/nCx, label = 'nC6x', color = '#e377c2', linestyle = ':')
plt.plot(T, nCm/nC, label = 'nCm', color = '#7f7f7f')
plt.plot(TEvolx, nCmx/nCx, label = 'nCmx', color = '#7f7f7f', linestyle = ':')
plt.yscale('log')
plt.xscale('log')
plt.ylim(2e-3, 1.2)
plt.xlim(1e4, 1e7)
plt.savefig("nC_vs_time.png", dpi = 300)
plt.close()

plt.plot(T, nN0/nN, label = 'nN0', color = '#1f77b4')
plt.plot(TEvolx, nN0x/nNx, label = 'nN0x', color = '#1f77b4', linestyle = ':')
plt.plot(T, nN1/nN, label = 'nN1', color = '#ff7f0e')
plt.plot(TEvolx, nN1x/nNx, label = 'nN1x', color = '#ff7f0e', linestyle = ':')
plt.plot(T, nN2/nN, label = 'nN2', color = '#2ca02c')
plt.plot(TEvolx, nN2x/nNx, label = 'nN2x', color = '#2ca02c', linestyle = ':')
plt.plot(T, nN3/nN, label = 'nN3', color = '#d62728')
plt.plot(TEvolx, nN3x/nNx, label = 'nN3x', color = '#d62728', linestyle = ':')
plt.plot(T, nN4/nN, label = 'nN4', color = '#9467bd')
plt.plot(TEvolx, nN4x/nNx, label = 'nN4x', color = '#9467bd', linestyle = ':')
plt.plot(T, nN5/nN, label = 'nN5', color = '#8c564b')
plt.plot(TEvolx, nN5x/nNx, label = 'nN5x', color = '#8c564b', linestyle = ':')
plt.plot(T, nN6/nN, label = 'nN6', color = '#e377c2')
plt.plot(TEvolx, nN6x/nNx, label = 'nN6x', color = '#e377c2', linestyle = ':')
plt.plot(T, nN7/nN, label = 'nN7', color = '#7f7f7f')
plt.plot(TEvolx, nN7x/nNx, label = 'nN7x', color = '#7f7f7f', linestyle = ':')
plt.yscale('log')
plt.xscale('log')
plt.ylim(2e-3, 1.2)
plt.xlim(1e4, 1e7)
plt.savefig("nN_vs_time.png", dpi = 300)
plt.close()

plt.plot(T, nO0/nO, label = 'nO0', color = '#1f77b4')
plt.plot(TEvolx, nO0x/nOx, label = 'nO0x', color = '#1f77b4', linestyle = ':')
plt.plot(T, nO1/nO, label = 'nO1', color = '#ff7f0e')
plt.plot(TEvolx, nO1x/nOx, label = 'nO1x', color = '#ff7f0e', linestyle = ':')
plt.plot(T, nO2/nO, label = 'nO2', color = '#2ca02c')
plt.plot(TEvolx, nO2x/nOx, label = 'nO2x', color = '#2ca02c', linestyle = ':')
plt.plot(T, nO3/nO, label = 'nO3', color = '#d62728')
plt.plot(TEvolx, nO3x/nOx, label = 'nO3x', color = '#d62728', linestyle = ':')
plt.plot(T, nO4/nO, label = 'nO4', color = '#9467bd')
plt.plot(TEvolx, nO4x/nOx, label = 'nO4x', color = '#9467bd', linestyle = ':')
plt.plot(T, nO5/nO, label = 'nO5', color = '#8c564b')
plt.plot(TEvolx, nO5x/nOx, label = 'nO5x', color = '#8c564b', linestyle = ':')
plt.plot(T, nO6/nO, label = 'nO6', color = '#e377c2')
plt.plot(TEvolx, nO6x/nOx, label = 'nO6x', color = '#e377c2', linestyle = ':')
plt.plot(T, nO7/nO, label = 'nO7', color = '#7f7f7f')
plt.plot(TEvolx, nO7x/nOx, label = 'nO7x', color = '#7f7f7f', linestyle = ':')
plt.plot(T, nO8/nO, label = 'nO8', color = '#bcbd22')
plt.plot(TEvolx, nO8x/nOx, label = 'nO8x', color = '#bcbd22', linestyle = ':')
plt.plot(T, nOm/nO, label = 'nOm', color = '#17becf')
plt.plot(TEvolx, nOmx/nOx, label = 'nOmx', color = '#17becf', linestyle = ':')
plt.yscale('log')
plt.xscale('log')
plt.ylim(2e-3, 1.2)
plt.xlim(1e4, 1e7)
plt.savefig("nO_vs_time.png", dpi = 300)
plt.close()

plt.plot(T, nNe0/nNe, label = 'nNe0', color = '#1f77b4')
plt.plot(TEvolx, nNe0x/nNex, label = 'nNe0x', color = '#1f77b4', linestyle = ':')
plt.plot(T, nNe1/nNe, label = 'nNe1', color = '#ff7f0e')
plt.plot(TEvolx, nNe1x/nNex, label = 'nNe1x', color = '#ff7f0e', linestyle = ':')
plt.plot(T, nNe2/nNe, label = 'nNe2', color = '#2ca02c')
plt.plot(TEvolx, nNe2x/nNex, label = 'nNe2x', color = '#2ca02c', linestyle = ':')
plt.plot(T, nNe3/nNe, label = 'nNe3', color = '#d62728')
plt.plot(TEvolx, nNe3x/nNex, label = 'nNe3x', color = '#d62728', linestyle = ':')
plt.plot(T, nNe4/nNe, label = 'nNe4', color = '#9467bd')
plt.plot(TEvolx, nNe4x/nNex, label = 'nNe4x', color = '#9467bd', linestyle = ':')
plt.plot(T, nNe5/nNe, label = 'nNe5', color = '#8c564b')
plt.plot(TEvolx, nNe5x/nNex, label = 'nNe5x', color = '#8c564b', linestyle = ':')
plt.plot(T, nNe6/nNe, label = 'nNe6', color = '#e377c2')
plt.plot(TEvolx, nNe6x/nNex, label = 'nNe6x', color = '#e377c2', linestyle = ':')
plt.plot(T, nNe7/nNe, label = 'nNe7', color = '#7f7f7f')
plt.plot(TEvolx, nNe7x/nNex, label = 'nNe7x', color = '#7f7f7f', linestyle = ':')
plt.plot(T, nNe8/nNe, label = 'nNe8', color = '#bcbd22')
plt.plot(TEvolx, nNe8x/nNex, label = 'nNe8x', color = '#bcbd22', linestyle = ':')
plt.plot(T, nNe9/nNe, label = 'nNe9', color = '#17becf')
plt.plot(TEvolx, nNe9x/nNex, label = 'nNe9x', color = '#17becf', linestyle = ':')
plt.plot(T, nNe10/nNe, label = 'nNe10', color = '#1f77b4')
plt.plot(TEvolx, nNe10x/nNex, label = 'nNe10x', color = '#1f77b4', linestyle = ':')
plt.yscale('log')
plt.xscale('log')
plt.ylim(2e-3, 1.2)
plt.xlim(1e4, 1e7)
plt.savefig("nNe_vs_time.png", dpi = 300)
plt.close()

plt.plot(T, nMg0/nMg, label = 'nMg0', color = '#1f77b4')
plt.plot(TEvolx, nMg0x/nMgx, label = 'nMg0x', color = '#1f77b4', linestyle = ':')
plt.plot(T, nMg1/nMg, label = 'nMg1', color = '#ff7f0e')
plt.plot(TEvolx, nMg1x/nMgx, label = 'nMg1x', color = '#ff7f0e', linestyle = ':')
plt.plot(T, nMg2/nMg, label = 'nMg2', color = '#2ca02c')
plt.plot(TEvolx, nMg2x/nMgx, label = 'nMg2x', color = '#2ca02c', linestyle = ':')
plt.plot(T, nMg3/nMg, label = 'nMg3', color = '#d62728')
plt.plot(TEvolx, nMg3x/nMgx, label = 'nMg3x', color = '#d62728', linestyle = ':')
plt.plot(T, nMg4/nMg, label = 'nMg4', color = '#9467bd')
plt.plot(TEvolx, nMg4x/nMgx, label = 'nMg4x', color = '#9467bd', linestyle = ':')
plt.plot(T, nMg5/nMg, label = 'nMg5', color = '#8c564b')
plt.plot(TEvolx, nMg5x/nMgx, label = 'nMg5x', color = '#8c564b', linestyle = ':')
plt.plot(T, nMg6/nMg, label = 'nMg6', color = '#e377c2')
plt.plot(TEvolx, nMg6x/nMgx, label = 'nMg6x', color = '#e377c2', linestyle = ':')
plt.plot(T, nMg7/nMg, label = 'nMg7', color = '#7f7f7f')
plt.plot(TEvolx, nMg7x/nMgx, label = 'nMg7x', color = '#7f7f7f', linestyle = ':')
plt.plot(T, nMg8/nMg, label = 'nMg8', color = '#bcbd22')
plt.plot(TEvolx, nMg8x/nMgx, label = 'nMg8x', color = '#bcbd22', linestyle = ':')
plt.plot(T, nMg9/nMg, label = 'nMg9', color = '#17becf')
plt.plot(TEvolx, nMg9x/nMgx, label = 'nMg9x', color = '#17becf', linestyle = ':')
plt.plot(T, nMg10/nMg, label = 'nMg10', color = '#1f77b4')
plt.plot(TEvolx, nMg10x/nMgx, label = 'nMg10x', color = '#1f77b4', linestyle = ':')
plt.plot(T, nMg11/nMg, label = 'nMg11', color = '#aec7e8')
plt.plot(TEvolx, nMg11x/nMgx, label = 'nMg11x', color = '#aec7e8', linestyle = ':')
plt.plot(T, nMg12/nMg, label = 'nMg12', color = '#ffbb78')
plt.plot(TEvolx, nMg12x/nMgx, label = 'nMg12x', color = '#ffbb78', linestyle = ':')
plt.yscale('log')
plt.xscale('log')
plt.ylim(2e-3, 1.2)
plt.xlim(1e4, 1e7)
plt.savefig("nMg_vs_time.png", dpi = 300)
plt.close()

plt.plot(T, nSi0/nSi, label = 'nSi0', color = '#1f77b4')
plt.plot(TEvolx, nSi0x/nSix, label = 'nSi0x', color = '#1f77b4', linestyle = ':')
plt.plot(T, nSi1/nSi, label = 'nSi1', color = '#ff7f0e')
plt.plot(TEvolx, nSi1x/nSix, label = 'nSi1x', color = '#ff7f0e', linestyle = ':')
plt.plot(T, nSi2/nSi, label = 'nSi2', color = '#2ca02c')
plt.plot(TEvolx, nSi2x/nSix, label = 'nSi2x', color = '#2ca02c', linestyle = ':')
plt.plot(T, nSi3/nSi, label = 'nSi3', color = '#d62728')
plt.plot(TEvolx, nSi3x/nSix, label = 'nSi3x', color = '#d62728', linestyle = ':')
plt.plot(T, nSi4/nSi, label = 'nSi4', color = '#9467bd')
plt.plot(TEvolx, nSi4x/nSix, label = 'nSi4x', color = '#9467bd', linestyle = ':')
plt.plot(T, nSi5/nSi, label = 'nSi5', color = '#8c564b')
plt.plot(TEvolx, nSi5x/nSix, label = 'nSi5x', color = '#8c564b', linestyle = ':')
plt.plot(T, nSi6/nSi, label = 'nSi6', color = '#e377c2')
plt.plot(TEvolx, nSi6x/nSix, label = 'nSi6x', color = '#e377c2', linestyle = ':')
plt.plot(T, nSi7/nSi, label = 'nSi7', color = '#7f7f7f')
plt.plot(TEvolx, nSi7x/nSix, label = 'nSi7x', color = '#7f7f7f', linestyle = ':')
plt.plot(T, nSi8/nSi, label = 'nSi8', color = '#bcbd22')
plt.plot(TEvolx, nSi8x/nSix, label = 'nSi8x', color = '#bcbd22', linestyle = ':')
plt.plot(T, nSi9/nSi, label = 'nSi9', color = '#17becf')
plt.plot(TEvolx, nSi9x/nSix, label = 'nSi9x', color = '#17becf', linestyle = ':')
plt.plot(T, nSi10/nSi, label = 'nSi10', color = '#1f77b4')
plt.plot(TEvolx, nSi10x/nSix, label = 'nSi10x', color = '#1f77b4', linestyle = ':')
plt.plot(T, nSi11/nSi, label = 'nSi11', color = '#aec7e8')
plt.plot(TEvolx, nSi11x/nSix, label = 'nSi11x', color = '#aec7e8', linestyle = ':')
plt.plot(T, nSi12/nSi, label = 'nSi12', color = '#ffbb78')
plt.plot(TEvolx, nSi12x/nSix, label = 'nSi12x', color = '#ffbb78', linestyle = ':')
plt.plot(T, nSi13/nSi, label = 'nSi13', color = '#98df8a')
plt.plot(TEvolx, nSi13x/nSix, label = 'nSi13x', color = '#98df8a', linestyle = ':')
plt.plot(T, nSi14/nSi, label = 'nSi14', color = '#ff9896')
plt.plot(TEvolx, nSi14x/nSix, label = 'nSi14x', color = '#ff9896', linestyle = ':')
plt.yscale('log')
plt.xscale('log')
plt.ylim(2e-3, 1.2)
plt.xlim(1e4, 1e7)
plt.savefig("nSi_vs_time.png", dpi = 300)
plt.close()

plt.plot(T, nS0/nS, label = 'nS0', color = '#1f77b4')
plt.plot(TEvolx, nS0x/nSx, label = 'nS0x', color = '#1f77b4', linestyle = ':')
plt.plot(T, nS1/nS, label = 'nS1', color = '#ff7f0e')
plt.plot(TEvolx, nS1x/nSx, label = 'nS1x', color = '#ff7f0e', linestyle = ':')
plt.plot(T, nS2/nS, label = 'nS2', color = '#2ca02c')
plt.plot(TEvolx, nS2x/nSx, label = 'nS2x', color = '#2ca02c', linestyle = ':')
plt.plot(T, nS3/nS, label = 'nS3', color = '#d62728')
plt.plot(TEvolx, nS3x/nSx, label = 'nS3x', color = '#d62728', linestyle = ':')
plt.plot(T, nS4/nS, label = 'nS4', color = '#9467bd')
plt.plot(TEvolx, nS4x/nSx, label = 'nS4x', color = '#9467bd', linestyle = ':')
plt.plot(T, nS5/nS, label = 'nS5', color = '#8c564b')
plt.plot(TEvolx, nS5x/nSx, label = 'nS5x', color = '#8c564b', linestyle = ':')
plt.plot(T, nS6/nS, label = 'nS6', color = '#e377c2')
plt.plot(TEvolx, nS6x/nSx, label = 'nS6x', color = '#e377c2', linestyle = ':')
plt.plot(T, nS7/nS, label = 'nS7', color = '#7f7f7f')
plt.plot(TEvolx, nS7x/nSx, label = 'nS7x', color = '#7f7f7f', linestyle = ':')
plt.plot(T, nS8/nS, label = 'nS8', color = '#bcbd22')
plt.plot(TEvolx, nS8x/nSx, label = 'nS8x', color = '#bcbd22', linestyle = ':')
plt.plot(T, nS9/nS, label = 'nS9', color = '#17becf')
plt.plot(TEvolx, nS9x/nSx, label = 'nS9x', color = '#17becf', linestyle = ':')
plt.plot(T, nS10/nS, label = 'nS10', color = '#1f77b4')
plt.plot(TEvolx, nS10x/nSx, label = 'nS10x', color = '#1f77b4', linestyle = ':')
plt.plot(T, nS11/nS, label = 'nS11', color = '#aec7e8')
plt.plot(TEvolx, nS11x/nSx, label = 'nS11x', color = '#aec7e8', linestyle = ':')
plt.plot(T, nS12/nS, label = 'nS12', color = '#ffbb78')
plt.plot(TEvolx, nS12x/nSx, label = 'nS12x', color = '#ffbb78', linestyle = ':')
plt.plot(T, nS13/nS, label = 'nS13', color = '#98df8a')
plt.plot(TEvolx, nS13x/nSx, label = 'nS13x', color = '#98df8a', linestyle = ':')
plt.plot(T, nS14/nS, label = 'nS14', color = '#ff9896')
plt.plot(TEvolx, nS14x/nSx, label = 'nS14x', color = '#ff9896', linestyle = ':')
plt.plot(T, nS15/nS, label = 'nS15', color = '#c5b0d5')
plt.plot(TEvolx, nS15x/nSx, label = 'nS15x', color = '#c5b0d5', linestyle = ':')
plt.plot(T, nS16/nS, label = 'nS16', color = '#c49c94')
plt.plot(TEvolx, nS16x/nSx, label = 'nS16x', color = '#c49c94', linestyle = ':')
plt.yscale('log')
plt.xscale('log')
plt.ylim(2e-3, 1.2)
plt.xlim(1e4, 1e7)
plt.savefig("nS_vs_time.png", dpi = 300)
plt.close()

print('Done !!!')

print('Elapsed time = ', time.time() - TA)






