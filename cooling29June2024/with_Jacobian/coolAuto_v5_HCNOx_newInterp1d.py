
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
           nOm):

  Tx = np.log10(T)

  ne = (
         1 * nHII - nHm + (nHeII + 2.0 * nHeIII) + 1 * nCII + 2 * nCIII + 3 * nCIV
       + 4 * nCV + 5 * nCVI + 6 * nCVII + -1 * nCm + 1 * nNII
       + 2 * nNIII + 3 * nNIV + 4 * nNV + 5 * nNVI + 6 * nNVII + 7 * nNVIII
       + 1 * nOII + 2 * nOIII + 3 * nOIV + 4 * nOV + 5 * nOVI
       + 6 * nOVII + 7 * nOVIII + 8 * nOIX + -1 * nOm
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
        + LCompton)

  return Lamb


#----- find_index
def find_index(Temp, T0):

  N_T = len(Temp)

  ndx_T = -1

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


#----- find_index
def interp1d_h(T, R, T0, ndxL):

  R0 = R[ndxL] + (R[ndxL + 1] - R[ndxL]) * (T0 - T[ndxL]) / (T[ndxL + 1] - T[ndxL])
  
  return R0


#----- func
def func(t, y):

  nHI, nHII, nHm, nHeI, nHeII, nHeIII, nCI, nCII, nCIII, nCIV, \
  nCV, nCVI, nCVII, nCm, nNI, nNII, nNIII, nNIV, nNV, nNVI, \
  nNVII, nNVIII, nOI, nOII, nOIII, nOIV, nOV, nOVI, nOVII, nOVIII, \
  nOIX, nOm, T = y

  Tx = np.log10(T)
  ndxL = find_index(Temp, Tx)
  
  
  R_HI_to_HII_via_e = 10**interp1d_h(Temp, R_HI_to_HII_via_e_, Tx, ndxL)
  R_HI_to_Hm_via_e = 10**interp1d_h(Temp, R_HI_to_Hm_via_e_, Tx, ndxL)
  R_Hm_to_HI_via_HI = 10**interp1d_h(Temp, R_Hm_to_HI_via_HI_, Tx, ndxL)
  R_HeII_to_HeI_via_HI = 10**interp1d_h(Temp, R_HeII_to_HeI_via_HI_, Tx, ndxL)
  R_HeIII_to_HeII_via_HI = 10**interp1d_h(Temp, R_HeIII_to_HeII_via_HI_, Tx, ndxL)
  R_CII_to_CI_via_HI = 10**interp1d_h(Temp, R_CII_to_CI_via_HI_, Tx, ndxL)
  R_CIII_to_CII_via_HI = 10**interp1d_h(Temp, R_CIII_to_CII_via_HI_, Tx, ndxL)
  R_CIV_to_CIII_via_HI = 10**interp1d_h(Temp, R_CIV_to_CIII_via_HI_, Tx, ndxL)
  R_CV_to_CIV_via_HI = 10**interp1d_h(Temp, R_CV_to_CIV_via_HI_, Tx, ndxL)
  R_CVI_to_CV_via_HI = 10**interp1d_h(Temp, R_CVI_to_CV_via_HI_, Tx, ndxL)
  R_NII_to_NI_via_HI = 10**interp1d_h(Temp, R_NII_to_NI_via_HI_, Tx, ndxL)
  R_NIII_to_NII_via_HI = 10**interp1d_h(Temp, R_NIII_to_NII_via_HI_, Tx, ndxL)
  R_NIV_to_NIII_via_HI = 10**interp1d_h(Temp, R_NIV_to_NIII_via_HI_, Tx, ndxL)
  R_NV_to_NIV_via_HI = 10**interp1d_h(Temp, R_NV_to_NIV_via_HI_, Tx, ndxL)
  R_NVI_to_NV_via_HI = 10**interp1d_h(Temp, R_NVI_to_NV_via_HI_, Tx, ndxL)
  R_OII_to_OI_via_HI = 10**interp1d_h(Temp, R_OII_to_OI_via_HI_, Tx, ndxL)
  R_OIII_to_OII_via_HI = 10**interp1d_h(Temp, R_OIII_to_OII_via_HI_, Tx, ndxL)
  R_OIV_to_OIII_via_HI = 10**interp1d_h(Temp, R_OIV_to_OIII_via_HI_, Tx, ndxL)
  R_OV_to_OIV_via_HI = 10**interp1d_h(Temp, R_OV_to_OIV_via_HI_, Tx, ndxL)
  R_OVI_to_OV_via_HI = 10**interp1d_h(Temp, R_OVI_to_OV_via_HI_, Tx, ndxL)
  R_Hm_to_HI_via_HII = 10**interp1d_h(Temp, R_Hm_to_HI_via_HII_, Tx, ndxL)
  R_HeI_to_HeII_via_HII = 10**interp1d_h(Temp, R_HeI_to_HeII_via_HII_, Tx, ndxL)
  R_CI_to_CII_via_HII = 10**interp1d_h(Temp, R_CI_to_CII_via_HII_, Tx, ndxL)
  R_Cm_to_CI_via_HII = 10**interp1d_h(Temp, R_Cm_to_CI_via_HII_, Tx, ndxL)
  R_NI_to_NII_via_HII = 10**interp1d_h(Temp, R_NI_to_NII_via_HII_, Tx, ndxL)
  R_OI_to_OII_via_HII = 10**interp1d_h(Temp, R_OI_to_OII_via_HII_, Tx, ndxL)
  R_Om_to_OI_via_HII = 10**interp1d_h(Temp, R_Om_to_OI_via_HII_, Tx, ndxL)
  R_Hm_to_HI_via_HI = 10**interp1d_h(Temp, R_Hm_to_HI_via_HI_, Tx, ndxL)
  R_Hm_to_HI_via_e = 10**interp1d_h(Temp, R_Hm_to_HI_via_e_, Tx, ndxL)
  R_Hm_to_HI_via_HII = 10**interp1d_h(Temp, R_Hm_to_HI_via_HII_, Tx, ndxL)
  R_HeII_to_HeI_via_Hm = 10**interp1d_h(Temp, R_HeII_to_HeI_via_Hm_, Tx, ndxL)
  R_HeI_to_HeII_via_HII = 10**interp1d_h(Temp, R_HeI_to_HeII_via_HII_, Tx, ndxL)
  R_HeI_to_HeII_via_e = 10**interp1d_h(Temp, R_HeI_to_HeII_via_e_, Tx, ndxL)
  R_CIV_to_CIII_via_HeI = 10**interp1d_h(Temp, R_CIV_to_CIII_via_HeI_, Tx, ndxL)
  R_CV_to_CIV_via_HeI = 10**interp1d_h(Temp, R_CV_to_CIV_via_HeI_, Tx, ndxL)
  R_NIII_to_NII_via_HeI = 10**interp1d_h(Temp, R_NIII_to_NII_via_HeI_, Tx, ndxL)
  R_NIV_to_NIII_via_HeI = 10**interp1d_h(Temp, R_NIV_to_NIII_via_HeI_, Tx, ndxL)
  R_NV_to_NIV_via_HeI = 10**interp1d_h(Temp, R_NV_to_NIV_via_HeI_, Tx, ndxL)
  R_OIII_to_OII_via_HeI = 10**interp1d_h(Temp, R_OIII_to_OII_via_HeI_, Tx, ndxL)
  R_OIV_to_OIII_via_HeI = 10**interp1d_h(Temp, R_OIV_to_OIII_via_HeI_, Tx, ndxL)
  R_OV_to_OIV_via_HeI = 10**interp1d_h(Temp, R_OV_to_OIV_via_HeI_, Tx, ndxL)
  R_HeII_to_HeIII_via_e = 10**interp1d_h(Temp, R_HeII_to_HeIII_via_e_, Tx, ndxL)
  R_HeII_to_HeI_via_Hm = 10**interp1d_h(Temp, R_HeII_to_HeI_via_Hm_, Tx, ndxL)
  R_HeII_to_HeI_via_HI = 10**interp1d_h(Temp, R_HeII_to_HeI_via_HI_, Tx, ndxL)
  R_CI_to_CII_via_HeII = 10**interp1d_h(Temp, R_CI_to_CII_via_HeII_, Tx, ndxL)
  R_CII_to_CIII_via_HeII = 10**interp1d_h(Temp, R_CII_to_CIII_via_HeII_, Tx, ndxL)
  R_NII_to_NIII_via_HeII = 10**interp1d_h(Temp, R_NII_to_NIII_via_HeII_, Tx, ndxL)
  R_OI_to_OII_via_HeII = 10**interp1d_h(Temp, R_OI_to_OII_via_HeII_, Tx, ndxL)
  R_HeIII_to_HeII_via_HI = 10**interp1d_h(Temp, R_HeIII_to_HeII_via_HI_, Tx, ndxL)
  R_HeIII_to_HeII_via_e = 10**interp1d_h(Temp, R_HeIII_to_HeII_via_e_, Tx, ndxL)
  R_CI_to_CII_via_HeII = 10**interp1d_h(Temp, R_CI_to_CII_via_HeII_, Tx, ndxL)
  R_CI_to_CII_via_HII = 10**interp1d_h(Temp, R_CI_to_CII_via_HII_, Tx, ndxL)
  R_CI_to_CII_via_e = 10**interp1d_h(Temp, R_CI_to_CII_via_e_, Tx, ndxL)
  R_CII_to_CI_via_HI = 10**interp1d_h(Temp, R_CII_to_CI_via_HI_, Tx, ndxL)
  R_CII_to_CIII_via_HeII = 10**interp1d_h(Temp, R_CII_to_CIII_via_HeII_, Tx, ndxL)
  R_CII_to_CI_via_e = 10**interp1d_h(Temp, R_CII_to_CI_via_e_, Tx, ndxL)
  R_CII_to_CIII_via_e = 10**interp1d_h(Temp, R_CII_to_CIII_via_e_, Tx, ndxL)
  R_CIII_to_CII_via_HI = 10**interp1d_h(Temp, R_CIII_to_CII_via_HI_, Tx, ndxL)
  R_CIII_to_CII_via_e = 10**interp1d_h(Temp, R_CIII_to_CII_via_e_, Tx, ndxL)
  R_CIII_to_CIV_via_e = 10**interp1d_h(Temp, R_CIII_to_CIV_via_e_, Tx, ndxL)
  R_CIV_to_CIII_via_HeI = 10**interp1d_h(Temp, R_CIV_to_CIII_via_HeI_, Tx, ndxL)
  R_CIV_to_CIII_via_HI = 10**interp1d_h(Temp, R_CIV_to_CIII_via_HI_, Tx, ndxL)
  R_CIV_to_CIII_via_e = 10**interp1d_h(Temp, R_CIV_to_CIII_via_e_, Tx, ndxL)
  R_CIV_to_CV_via_e = 10**interp1d_h(Temp, R_CIV_to_CV_via_e_, Tx, ndxL)
  R_CV_to_CVI_via_e = 10**interp1d_h(Temp, R_CV_to_CVI_via_e_, Tx, ndxL)
  R_CV_to_CIV_via_e = 10**interp1d_h(Temp, R_CV_to_CIV_via_e_, Tx, ndxL)
  R_CV_to_CIV_via_HI = 10**interp1d_h(Temp, R_CV_to_CIV_via_HI_, Tx, ndxL)
  R_CV_to_CIV_via_HeI = 10**interp1d_h(Temp, R_CV_to_CIV_via_HeI_, Tx, ndxL)
  R_CVI_to_CVII_via_e = 10**interp1d_h(Temp, R_CVI_to_CVII_via_e_, Tx, ndxL)
  R_CVI_to_CV_via_HI = 10**interp1d_h(Temp, R_CVI_to_CV_via_HI_, Tx, ndxL)
  R_CVI_to_CV_via_e = 10**interp1d_h(Temp, R_CVI_to_CV_via_e_, Tx, ndxL)
  R_CVII_to_CVI_via_e = 10**interp1d_h(Temp, R_CVII_to_CVI_via_e_, Tx, ndxL)
  R_Cm_to_CI_via_HII = 10**interp1d_h(Temp, R_Cm_to_CI_via_HII_, Tx, ndxL)
  R_NI_to_NII_via_HII = 10**interp1d_h(Temp, R_NI_to_NII_via_HII_, Tx, ndxL)
  R_NI_to_NII_via_e = 10**interp1d_h(Temp, R_NI_to_NII_via_e_, Tx, ndxL)
  R_NII_to_NI_via_HI = 10**interp1d_h(Temp, R_NII_to_NI_via_HI_, Tx, ndxL)
  R_NII_to_NIII_via_HeII = 10**interp1d_h(Temp, R_NII_to_NIII_via_HeII_, Tx, ndxL)
  R_NII_to_NI_via_e = 10**interp1d_h(Temp, R_NII_to_NI_via_e_, Tx, ndxL)
  R_NII_to_NIII_via_e = 10**interp1d_h(Temp, R_NII_to_NIII_via_e_, Tx, ndxL)
  R_NIII_to_NII_via_HeI = 10**interp1d_h(Temp, R_NIII_to_NII_via_HeI_, Tx, ndxL)
  R_NIII_to_NII_via_HI = 10**interp1d_h(Temp, R_NIII_to_NII_via_HI_, Tx, ndxL)
  R_NIII_to_NII_via_e = 10**interp1d_h(Temp, R_NIII_to_NII_via_e_, Tx, ndxL)
  R_NIII_to_NIV_via_e = 10**interp1d_h(Temp, R_NIII_to_NIV_via_e_, Tx, ndxL)
  R_NIV_to_NIII_via_HeI = 10**interp1d_h(Temp, R_NIV_to_NIII_via_HeI_, Tx, ndxL)
  R_NIV_to_NIII_via_HI = 10**interp1d_h(Temp, R_NIV_to_NIII_via_HI_, Tx, ndxL)
  R_NIV_to_NIII_via_e = 10**interp1d_h(Temp, R_NIV_to_NIII_via_e_, Tx, ndxL)
  R_NIV_to_NV_via_e = 10**interp1d_h(Temp, R_NIV_to_NV_via_e_, Tx, ndxL)
  R_NV_to_NIV_via_HeI = 10**interp1d_h(Temp, R_NV_to_NIV_via_HeI_, Tx, ndxL)
  R_NV_to_NIV_via_HI = 10**interp1d_h(Temp, R_NV_to_NIV_via_HI_, Tx, ndxL)
  R_NV_to_NIV_via_e = 10**interp1d_h(Temp, R_NV_to_NIV_via_e_, Tx, ndxL)
  R_NV_to_NVI_via_e = 10**interp1d_h(Temp, R_NV_to_NVI_via_e_, Tx, ndxL)
  R_NVI_to_NV_via_HI = 10**interp1d_h(Temp, R_NVI_to_NV_via_HI_, Tx, ndxL)
  R_NVI_to_NVII_via_e = 10**interp1d_h(Temp, R_NVI_to_NVII_via_e_, Tx, ndxL)
  R_NVI_to_NV_via_e = 10**interp1d_h(Temp, R_NVI_to_NV_via_e_, Tx, ndxL)
  R_NVII_to_NVI_via_e = 10**interp1d_h(Temp, R_NVII_to_NVI_via_e_, Tx, ndxL)
  R_NVII_to_NVIII_via_e = 10**interp1d_h(Temp, R_NVII_to_NVIII_via_e_, Tx, ndxL)
  R_NVIII_to_NVII_via_e = 10**interp1d_h(Temp, R_NVIII_to_NVII_via_e_, Tx, ndxL)
  R_OI_to_OII_via_HeII = 10**interp1d_h(Temp, R_OI_to_OII_via_HeII_, Tx, ndxL)
  R_OI_to_OII_via_e = 10**interp1d_h(Temp, R_OI_to_OII_via_e_, Tx, ndxL)
  R_OI_to_OII_via_HII = 10**interp1d_h(Temp, R_OI_to_OII_via_HII_, Tx, ndxL)
  R_OII_to_OI_via_HI = 10**interp1d_h(Temp, R_OII_to_OI_via_HI_, Tx, ndxL)
  R_OII_to_OI_via_e = 10**interp1d_h(Temp, R_OII_to_OI_via_e_, Tx, ndxL)
  R_OII_to_OIII_via_e = 10**interp1d_h(Temp, R_OII_to_OIII_via_e_, Tx, ndxL)
  R_OIII_to_OII_via_HeI = 10**interp1d_h(Temp, R_OIII_to_OII_via_HeI_, Tx, ndxL)
  R_OIII_to_OII_via_HI = 10**interp1d_h(Temp, R_OIII_to_OII_via_HI_, Tx, ndxL)
  R_OIII_to_OII_via_e = 10**interp1d_h(Temp, R_OIII_to_OII_via_e_, Tx, ndxL)
  R_OIII_to_OIV_via_e = 10**interp1d_h(Temp, R_OIII_to_OIV_via_e_, Tx, ndxL)
  R_OIV_to_OV_via_e = 10**interp1d_h(Temp, R_OIV_to_OV_via_e_, Tx, ndxL)
  R_OIV_to_OIII_via_e = 10**interp1d_h(Temp, R_OIV_to_OIII_via_e_, Tx, ndxL)
  R_OIV_to_OIII_via_HI = 10**interp1d_h(Temp, R_OIV_to_OIII_via_HI_, Tx, ndxL)
  R_OIV_to_OIII_via_HeI = 10**interp1d_h(Temp, R_OIV_to_OIII_via_HeI_, Tx, ndxL)
  R_OV_to_OIV_via_HeI = 10**interp1d_h(Temp, R_OV_to_OIV_via_HeI_, Tx, ndxL)
  R_OV_to_OIV_via_HI = 10**interp1d_h(Temp, R_OV_to_OIV_via_HI_, Tx, ndxL)
  R_OV_to_OIV_via_e = 10**interp1d_h(Temp, R_OV_to_OIV_via_e_, Tx, ndxL)
  R_OV_to_OVI_via_e = 10**interp1d_h(Temp, R_OV_to_OVI_via_e_, Tx, ndxL)
  R_OVI_to_OV_via_HI = 10**interp1d_h(Temp, R_OVI_to_OV_via_HI_, Tx, ndxL)
  R_OVI_to_OV_via_e = 10**interp1d_h(Temp, R_OVI_to_OV_via_e_, Tx, ndxL)
  R_OVI_to_OVII_via_e = 10**interp1d_h(Temp, R_OVI_to_OVII_via_e_, Tx, ndxL)
  R_OVII_to_OVI_via_e = 10**interp1d_h(Temp, R_OVII_to_OVI_via_e_, Tx, ndxL)
  R_OVII_to_OVIII_via_e = 10**interp1d_h(Temp, R_OVII_to_OVIII_via_e_, Tx, ndxL)
  R_OVIII_to_OVII_via_e = 10**interp1d_h(Temp, R_OVIII_to_OVII_via_e_, Tx, ndxL)
  R_OVIII_to_OIX_via_e = 10**interp1d_h(Temp, R_OVIII_to_OIX_via_e_, Tx, ndxL)
  R_OIX_to_OVIII_via_e = 10**interp1d_h(Temp, R_OIX_to_OVIII_via_e_, Tx, ndxL)
  R_Om_to_OI_via_HII = 10**interp1d_h(Temp, R_Om_to_OI_via_HII_, Tx, ndxL)

  ne = (
       + 1 * nHII + -1 * nHm + 1 * nHeII + 2 * nHeIII + 1 * nCII + 2 * nCIII
       + 3 * nCIV + 4 * nCV + 5 * nCVI + 6 * nCVII + -1 * nCm + 1 * nNII
       + 2 * nNIII + 3 * nNIV + 4 * nNV + 5 * nNVI + 6 * nNVII + 7 * nNVIII
       + 1 * nOII + 2 * nOIII + 3 * nOIV + 4 * nOV + 5 * nOVI + 6 * nOVII
       + 7 * nOVIII + 8 * nOIX + -1 * nOm
       )

  ntot = (
           ne + nHI + nHII + nHm + nHeI + nHeII + nHeIII
         + nCI + nCII + nCIII + nCIV + nCV + nCVI
         + nCVII + nCm + nNI + nNII + nNIII + nNIV
         + nNV + nNVI + nNVII + nNVIII + nOI + nOII
         + nOIII + nOIV + nOV + nOVI + nOVII + nOVIII
         + nOIX + nOm
         )

  grain_rec_HII_to_HI = grain_recomb_rate("HII", T, ne, G0, A_v, Temp, Psi)
  grain_rec_HeII_to_HeI = grain_recomb_rate("HeII", T, ne, G0, A_v, Temp, Psi)
  grain_rec_CII_to_CI = grain_recomb_rate("CII", T, ne, G0, A_v, Temp, Psi)
  grain_rec_OII_to_OI = grain_recomb_rate("OII", T, ne, G0, A_v, Temp, Psi)

  const_OI_e_to_Om_ = 1.5000E-15 # constant/rates
  const_CI_e_to_Cm_ = 2.2500E-15 # constant/rates

  dnHI_dt = (
             - R_HI_to_HII_via_e * nHI * ne
             + 10**grain_rec_HII_to_HI * nHII * ne # grain_recombination
             - R_HI_to_Hm_via_e * nHI * ne
             - R_Hm_to_HI_via_HI * nHm * nHI
             + R_Hm_to_HI_via_HI * nHm * nHI
             + R_Hm_to_HI_via_e * nHm * ne
             + R_Hm_to_HI_via_HII * nHm * nHII
             + R_HeI_to_HeII_via_HII * nHeI * nHII
             + R_HeII_to_HeI_via_Hm * nHeII * nHm
             - R_HeII_to_HeI_via_HI * nHeII * nHI
             - R_HeIII_to_HeII_via_HI * nHeIII * nHI
             + R_CI_to_CII_via_HII * nCI * nHII
             - R_CII_to_CI_via_HI * nCII * nHI
             - R_CIII_to_CII_via_HI * nCIII * nHI
             - R_CIV_to_CIII_via_HI * nCIV * nHI
             - R_CV_to_CIV_via_HI * nCV * nHI
             - R_CVI_to_CV_via_HI * nCVI * nHI
             + R_Cm_to_CI_via_HII * nCm * nHII
             + R_NI_to_NII_via_HII * nNI * nHII
             - R_NII_to_NI_via_HI * nNII * nHI
             - R_NIII_to_NII_via_HI * nNIII * nHI
             - R_NIV_to_NIII_via_HI * nNIV * nHI
             - R_NV_to_NIV_via_HI * nNV * nHI
             - R_NVI_to_NV_via_HI * nNVI * nHI
             + R_OI_to_OII_via_HII * nOI * nHII
             - R_OII_to_OI_via_HI * nOII * nHI
             - R_OIII_to_OII_via_HI * nOIII * nHI
             - R_OIV_to_OIII_via_HI * nOIV * nHI
             - R_OV_to_OIV_via_HI * nOV * nHI
             - R_OVI_to_OV_via_HI * nOVI * nHI
             + R_Om_to_OI_via_HII * nOm * nHII
             + 10**R_HII_to_HI_via_e_caseA(Tx) * nHII * ne # H CaseA
            )

  dnHII_dt = (
              + R_HI_to_HII_via_e * nHI * ne
              - 10**grain_rec_HII_to_HI * nHII * ne # grain_recombination
              - R_Hm_to_HI_via_HII * nHm * nHII
              - R_HeI_to_HeII_via_HII * nHeI * nHII
              + R_HeII_to_HeI_via_HI * nHeII * nHI
              + R_HeIII_to_HeII_via_HI * nHeIII * nHI
              - R_CI_to_CII_via_HII * nCI * nHII
              + R_CII_to_CI_via_HI * nCII * nHI
              + R_CIII_to_CII_via_HI * nCIII * nHI
              + R_CIV_to_CIII_via_HI * nCIV * nHI
              + R_CV_to_CIV_via_HI * nCV * nHI
              + R_CVI_to_CV_via_HI * nCVI * nHI
              - R_Cm_to_CI_via_HII * nCm * nHII
              - R_NI_to_NII_via_HII * nNI * nHII
              + R_NII_to_NI_via_HI * nNII * nHI
              + R_NIII_to_NII_via_HI * nNIII * nHI
              + R_NIV_to_NIII_via_HI * nNIV * nHI
              + R_NV_to_NIV_via_HI * nNV * nHI
              + R_NVI_to_NV_via_HI * nNVI * nHI
              - R_OI_to_OII_via_HII * nOI * nHII
              + R_OII_to_OI_via_HI * nOII * nHI
              + R_OIII_to_OII_via_HI * nOIII * nHI
              + R_OIV_to_OIII_via_HI * nOIV * nHI
              + R_OV_to_OIV_via_HI * nOV * nHI
              + R_OVI_to_OV_via_HI * nOVI * nHI
              - R_Om_to_OI_via_HII * nOm * nHII
              - 10**R_HII_to_HI_via_e_caseA(Tx) * nHII * ne # H CaseA
             )

  dnHm_dt = (
             + R_HI_to_Hm_via_e * nHI * ne
             - R_Hm_to_HI_via_HI * nHm * nHI
             - R_Hm_to_HI_via_e * nHm * ne
             - R_Hm_to_HI_via_HII * nHm * nHII
             - R_HeII_to_HeI_via_Hm * nHeII * nHm
            )

  dnHeI_dt = (
              - R_HeI_to_HeII_via_HII * nHeI * nHII
              + 10**grain_rec_HeII_to_HeI * nHeII * ne # grain_recombination
              - R_HeI_to_HeII_via_e * nHeI * ne
              + R_HeII_to_HeI_via_Hm * nHeII * nHm
              + R_HeII_to_HeI_via_HI * nHeII * nHI
              + R_CI_to_CII_via_HeII * nCI * nHeII
              + R_CII_to_CIII_via_HeII * nCII * nHeII
              - R_CIV_to_CIII_via_HeI * nCIV * nHeI
              - R_CV_to_CIV_via_HeI * nCV * nHeI
              + R_NII_to_NIII_via_HeII * nNII * nHeII
              - R_NIII_to_NII_via_HeI * nNIII * nHeI
              - R_NIV_to_NIII_via_HeI * nNIV * nHeI
              - R_NV_to_NIV_via_HeI * nNV * nHeI
              + R_OI_to_OII_via_HeII * nOI * nHeII
              - R_OIII_to_OII_via_HeI * nOIII * nHeI
              - R_OIV_to_OIII_via_HeI * nOIV * nHeI
              - R_OV_to_OIV_via_HeI * nOV * nHeI
              + 10**R_HeII_to_HeI_via_e_caseA(Tx) * nHeII * ne # He CaseA
             )

  dnHeII_dt = (
               + R_HeI_to_HeII_via_HII * nHeI * nHII
               - 10**grain_rec_HeII_to_HeI * nHeII * ne # grain_recombination
               + R_HeI_to_HeII_via_e * nHeI * ne
               - R_HeII_to_HeIII_via_e * nHeII * ne
               - R_HeII_to_HeI_via_Hm * nHeII * nHm
               - R_HeII_to_HeI_via_HI * nHeII * nHI
               + R_HeIII_to_HeII_via_HI * nHeIII * nHI
               + R_HeIII_to_HeII_via_e * nHeIII * ne
               - R_CI_to_CII_via_HeII * nCI * nHeII
               - R_CII_to_CIII_via_HeII * nCII * nHeII
               + R_CIV_to_CIII_via_HeI * nCIV * nHeI
               + R_CV_to_CIV_via_HeI * nCV * nHeI
               - R_NII_to_NIII_via_HeII * nNII * nHeII
               + R_NIII_to_NII_via_HeI * nNIII * nHeI
               + R_NIV_to_NIII_via_HeI * nNIV * nHeI
               + R_NV_to_NIV_via_HeI * nNV * nHeI
               - R_OI_to_OII_via_HeII * nOI * nHeII
               + R_OIII_to_OII_via_HeI * nOIII * nHeI
               + R_OIV_to_OIII_via_HeI * nOIV * nHeI
               + R_OV_to_OIV_via_HeI * nOV * nHeI
               - 10**R_HeII_to_HeI_via_e_caseA(Tx) * nHeII * ne # He CaseA
              )

  dnHeIII_dt = (
                + R_HeII_to_HeIII_via_e * nHeII * ne
                - R_HeIII_to_HeII_via_HI * nHeIII * nHI
                - R_HeIII_to_HeII_via_e * nHeIII * ne
               )

  dnCI_dt = (
             - R_CI_to_CII_via_HeII * nCI * nHeII
             + 10**grain_rec_CII_to_CI * nCII * ne # grain_recombination
             - R_CI_to_CII_via_HII * nCI * nHII
             - R_CI_to_CII_via_e * nCI * ne
             + R_CII_to_CI_via_HI * nCII * nHI
             + R_CII_to_CI_via_e * nCII * ne
             + R_Cm_to_CI_via_HII * nCm * nHII
             - const_CI_e_to_Cm_ * nCI * ne # constant rate
            )

  dnCII_dt = (
              + R_CI_to_CII_via_HeII * nCI * nHeII
              - 10**grain_rec_CII_to_CI * nCII * ne # grain_recombination
              + R_CI_to_CII_via_HII * nCI * nHII
              + R_CI_to_CII_via_e * nCI * ne
              - R_CII_to_CI_via_HI * nCII * nHI
              - R_CII_to_CIII_via_HeII * nCII * nHeII
              - R_CII_to_CI_via_e * nCII * ne
              - R_CII_to_CIII_via_e * nCII * ne
              + R_CIII_to_CII_via_HI * nCIII * nHI
              + R_CIII_to_CII_via_e * nCIII * ne
             )

  dnCIII_dt = (
               + R_CII_to_CIII_via_HeII * nCII * nHeII
               + R_CII_to_CIII_via_e * nCII * ne
               - R_CIII_to_CII_via_HI * nCIII * nHI
               - R_CIII_to_CII_via_e * nCIII * ne
               - R_CIII_to_CIV_via_e * nCIII * ne
               + R_CIV_to_CIII_via_HeI * nCIV * nHeI
               + R_CIV_to_CIII_via_HI * nCIV * nHI
               + R_CIV_to_CIII_via_e * nCIV * ne
              )

  dnCIV_dt = (
              + R_CIII_to_CIV_via_e * nCIII * ne
              - R_CIV_to_CIII_via_HeI * nCIV * nHeI
              - R_CIV_to_CIII_via_HI * nCIV * nHI
              - R_CIV_to_CIII_via_e * nCIV * ne
              - R_CIV_to_CV_via_e * nCIV * ne
              + R_CV_to_CIV_via_e * nCV * ne
              + R_CV_to_CIV_via_HI * nCV * nHI
              + R_CV_to_CIV_via_HeI * nCV * nHeI
             )

  dnCV_dt = (
             + R_CIV_to_CV_via_e * nCIV * ne
             - R_CV_to_CVI_via_e * nCV * ne
             - R_CV_to_CIV_via_e * nCV * ne
             - R_CV_to_CIV_via_HI * nCV * nHI
             - R_CV_to_CIV_via_HeI * nCV * nHeI
             + R_CVI_to_CV_via_HI * nCVI * nHI
             + R_CVI_to_CV_via_e * nCVI * ne
            )

  dnCVI_dt = (
              + R_CV_to_CVI_via_e * nCV * ne
              - R_CVI_to_CVII_via_e * nCVI * ne
              - R_CVI_to_CV_via_HI * nCVI * nHI
              - R_CVI_to_CV_via_e * nCVI * ne
              + R_CVII_to_CVI_via_e * nCVII * ne
             )

  dnCVII_dt = (
               + R_CVI_to_CVII_via_e * nCVI * ne
               - R_CVII_to_CVI_via_e * nCVII * ne
              )

  dnCm_dt = (
             - R_Cm_to_CI_via_HII * nCm * nHII
             + const_CI_e_to_Cm_ * nCI * ne # constant rate
            )

  dnNI_dt = (
             - R_NI_to_NII_via_HII * nNI * nHII
             - R_NI_to_NII_via_e * nNI * ne
             + R_NII_to_NI_via_HI * nNII * nHI
             + R_NII_to_NI_via_e * nNII * ne
            )

  dnNII_dt = (
              + R_NI_to_NII_via_HII * nNI * nHII
              + R_NI_to_NII_via_e * nNI * ne
              - R_NII_to_NI_via_HI * nNII * nHI
              - R_NII_to_NIII_via_HeII * nNII * nHeII
              - R_NII_to_NI_via_e * nNII * ne
              - R_NII_to_NIII_via_e * nNII * ne
              + R_NIII_to_NII_via_HeI * nNIII * nHeI
              + R_NIII_to_NII_via_HI * nNIII * nHI
              + R_NIII_to_NII_via_e * nNIII * ne
             )

  dnNIII_dt = (
               + R_NII_to_NIII_via_HeII * nNII * nHeII
               + R_NII_to_NIII_via_e * nNII * ne
               - R_NIII_to_NII_via_HeI * nNIII * nHeI
               - R_NIII_to_NII_via_HI * nNIII * nHI
               - R_NIII_to_NII_via_e * nNIII * ne
               - R_NIII_to_NIV_via_e * nNIII * ne
               + R_NIV_to_NIII_via_HeI * nNIV * nHeI
               + R_NIV_to_NIII_via_HI * nNIV * nHI
               + R_NIV_to_NIII_via_e * nNIV * ne
              )

  dnNIV_dt = (
              + R_NIII_to_NIV_via_e * nNIII * ne
              - R_NIV_to_NIII_via_HeI * nNIV * nHeI
              - R_NIV_to_NIII_via_HI * nNIV * nHI
              - R_NIV_to_NIII_via_e * nNIV * ne
              - R_NIV_to_NV_via_e * nNIV * ne
              + R_NV_to_NIV_via_HeI * nNV * nHeI
              + R_NV_to_NIV_via_HI * nNV * nHI
              + R_NV_to_NIV_via_e * nNV * ne
             )

  dnNV_dt = (
             + R_NIV_to_NV_via_e * nNIV * ne
             - R_NV_to_NIV_via_HeI * nNV * nHeI
             - R_NV_to_NIV_via_HI * nNV * nHI
             - R_NV_to_NIV_via_e * nNV * ne
             - R_NV_to_NVI_via_e * nNV * ne
             + R_NVI_to_NV_via_HI * nNVI * nHI
             + R_NVI_to_NV_via_e * nNVI * ne
            )

  dnNVI_dt = (
              + R_NV_to_NVI_via_e * nNV * ne
              - R_NVI_to_NV_via_HI * nNVI * nHI
              - R_NVI_to_NVII_via_e * nNVI * ne
              - R_NVI_to_NV_via_e * nNVI * ne
              + R_NVII_to_NVI_via_e * nNVII * ne
             )

  dnNVII_dt = (
               + R_NVI_to_NVII_via_e * nNVI * ne
               - R_NVII_to_NVI_via_e * nNVII * ne
               - R_NVII_to_NVIII_via_e * nNVII * ne
               + R_NVIII_to_NVII_via_e * nNVIII * ne
              )

  dnNVIII_dt = (
                + R_NVII_to_NVIII_via_e * nNVII * ne
                - R_NVIII_to_NVII_via_e * nNVIII * ne
               )

  dnOI_dt = (
             - R_OI_to_OII_via_HeII * nOI * nHeII
             + 10**grain_rec_OII_to_OI * nOII * ne # grain_recombination
             - R_OI_to_OII_via_e * nOI * ne
             - R_OI_to_OII_via_HII * nOI * nHII
             + R_OII_to_OI_via_HI * nOII * nHI
             + R_OII_to_OI_via_e * nOII * ne
             + R_Om_to_OI_via_HII * nOm * nHII
             - const_OI_e_to_Om_ * nOI * ne # constant rate
            )

  dnOII_dt = (
              + R_OI_to_OII_via_HeII * nOI * nHeII
              - 10**grain_rec_OII_to_OI * nOII * ne # grain_recombination
              + R_OI_to_OII_via_e * nOI * ne
              + R_OI_to_OII_via_HII * nOI * nHII
              - R_OII_to_OI_via_HI * nOII * nHI
              - R_OII_to_OI_via_e * nOII * ne
              - R_OII_to_OIII_via_e * nOII * ne
              + R_OIII_to_OII_via_HeI * nOIII * nHeI
              + R_OIII_to_OII_via_HI * nOIII * nHI
              + R_OIII_to_OII_via_e * nOIII * ne
             )

  dnOIII_dt = (
               + R_OII_to_OIII_via_e * nOII * ne
               - R_OIII_to_OII_via_HeI * nOIII * nHeI
               - R_OIII_to_OII_via_HI * nOIII * nHI
               - R_OIII_to_OII_via_e * nOIII * ne
               - R_OIII_to_OIV_via_e * nOIII * ne
               + R_OIV_to_OIII_via_e * nOIV * ne
               + R_OIV_to_OIII_via_HI * nOIV * nHI
               + R_OIV_to_OIII_via_HeI * nOIV * nHeI
              )

  dnOIV_dt = (
              + R_OIII_to_OIV_via_e * nOIII * ne
              - R_OIV_to_OV_via_e * nOIV * ne
              - R_OIV_to_OIII_via_e * nOIV * ne
              - R_OIV_to_OIII_via_HI * nOIV * nHI
              - R_OIV_to_OIII_via_HeI * nOIV * nHeI
              + R_OV_to_OIV_via_HeI * nOV * nHeI
              + R_OV_to_OIV_via_HI * nOV * nHI
              + R_OV_to_OIV_via_e * nOV * ne
             )

  dnOV_dt = (
             + R_OIV_to_OV_via_e * nOIV * ne
             - R_OV_to_OIV_via_HeI * nOV * nHeI
             - R_OV_to_OIV_via_HI * nOV * nHI
             - R_OV_to_OIV_via_e * nOV * ne
             - R_OV_to_OVI_via_e * nOV * ne
             + R_OVI_to_OV_via_HI * nOVI * nHI
             + R_OVI_to_OV_via_e * nOVI * ne
            )

  dnOVI_dt = (
              + R_OV_to_OVI_via_e * nOV * ne
              - R_OVI_to_OV_via_HI * nOVI * nHI
              - R_OVI_to_OV_via_e * nOVI * ne
              - R_OVI_to_OVII_via_e * nOVI * ne
              + R_OVII_to_OVI_via_e * nOVII * ne
             )

  dnOVII_dt = (
               + R_OVI_to_OVII_via_e * nOVI * ne
               - R_OVII_to_OVI_via_e * nOVII * ne
               - R_OVII_to_OVIII_via_e * nOVII * ne
               + R_OVIII_to_OVII_via_e * nOVIII * ne
              )

  dnOVIII_dt = (
                + R_OVII_to_OVIII_via_e * nOVII * ne
                - R_OVIII_to_OVII_via_e * nOVIII * ne
                - R_OVIII_to_OIX_via_e * nOVIII * ne
                + R_OIX_to_OVIII_via_e * nOIX * ne
               )

  dnOIX_dt = (
              + R_OVIII_to_OIX_via_e * nOVIII * ne
              - R_OIX_to_OVIII_via_e * nOIX * ne
             )

  dnOm_dt = (
             - R_Om_to_OI_via_HII * nOm * nHII
             + const_OI_e_to_Om_ * nOI * ne # constant rate
            )

  Lamb = Lambda(
                T, nHI, nHII, nHm, nHeI, nHeII, nHeIII, nCI, nCII, nCIII, nCIV, 
                nCV, nCVI, nCVII, nCm, nNI, nNII, nNIII, nNIV, nNV, nNVI, 
                nNVII, nNVIII, nOI, nOII, nOIII, nOIV, nOV, nOVI, nOVII, nOVIII, 
                nOIX, nOm)

  dne_dt = (
           + 1 * dnHII_dt + -1 * dnHm_dt + 1 * dnHeII_dt + 2 * dnHeIII_dt + 1 * dnCII_dt + 2 * dnCIII_dt
           + 3 * dnCIV_dt + 4 * dnCV_dt + 5 * dnCVI_dt + 6 * dnCVII_dt + -1 * dnCm_dt + 1 * dnNII_dt
           + 2 * dnNIII_dt + 3 * dnNIV_dt + 4 * dnNV_dt + 5 * dnNVI_dt + 6 * dnNVII_dt + 7 * dnNVIII_dt
           + 1 * dnOII_dt + 2 * dnOIII_dt + 3 * dnOIV_dt + 4 * dnOV_dt + 5 * dnOVI_dt + 6 * dnOVII_dt
           + 7 * dnOVIII_dt + 8 * dnOIX_dt + -1 * dnOm_dt
           )

  dntot_dt = (
                dne_dt + dnHI_dt + dnHII_dt + dnHm_dt + dnHeI_dt + dnHeII_dt + dnHeIII_dt
              + dnCI_dt + dnCII_dt + dnCIII_dt + dnCIV_dt + dnCV_dt + dnCVI_dt
              + dnCVII_dt + dnCm_dt + dnNI_dt + dnNII_dt + dnNIII_dt + dnNIV_dt
              + dnNV_dt + dnNVI_dt + dnNVII_dt + dnNVIII_dt + dnOI_dt + dnOII_dt
              + dnOIII_dt + dnOIV_dt + dnOV_dt + dnOVI_dt + dnOVII_dt + dnOVIII_dt
              + dnOIX_dt + dnOm_dt
             )

  dT_dt = -1.0 * (gamma - 1.0) / kB / ntot * (Lamb + 1. / (gamma - 1.) * kB * T * dntot_dt)

  return [
          dnHI_dt, dnHII_dt, dnHm_dt, dnHeI_dt, dnHeII_dt, dnHeIII_dt,
          dnCI_dt, dnCII_dt, dnCIII_dt, dnCIV_dt, dnCV_dt, dnCVI_dt,
          dnCVII_dt, dnCm_dt, dnNI_dt, dnNII_dt, dnNIII_dt, dnNIV_dt,
          dnNV_dt, dnNVI_dt, dnNVII_dt, dnNVIII_dt, dnOI_dt, dnOII_dt,
          dnOIII_dt, dnOIV_dt, dnOV_dt, dnOVI_dt, dnOVII_dt, dnOVIII_dt,
          dnOIX_dt, dnOm_dt, dT_dt
         ]







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

y0 = [
      nH0_i, nHp_i, nHm_i, nHe0_i, nHep_i, nHepp_i, 
      nC0_i, nC1_i, nC2_i, nC3_i, nC4_i, nC5_i, 
      nC6_i, nCm_i, nN0_i, nN1_i, nN2_i, nN3_i, 
      nN4_i, nN5_i, nN6_i, nN7_i, nO0_i, nO1_i, 
      nO2_i, nO3_i, nO4_i, nO5_i, nO6_i, nO7_i, 
      nO8_i, nOm_i, 
      T_i
      ]

A_v = 1.0
G0 = 0.01
dust_ratio = 0.01


t_span = (1*3.16e7, 10000*3.16e7)

T001 = time.time()
solution = solve_ivp(func, t_span, y0, method="LSODA", dense_output=True)
print(f'Elapsed time = {time.time() - T001}')
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


T = y[32, :]

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

print('Done !!!')



