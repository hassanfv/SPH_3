
import numpy as np
import matplotlib.pyplot as plt
#from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp
import pickle
from datax_vz_t import *
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
           nSXI, nSXII, nSXIII, nSXIV, nSXV, nSXVI, nSXVII, nCaI, nCaII, nCaIII, 
           nCaIV, nCaV, nCaVI, nCaVII, nCaVIII, nCaIX, nCaX, nCaXI, nCaXII, nCaXIII, 
           nCaXIV, nCaXV, nCaXVI, nCaXVII, nCaXVIII, nCaXIX, nCaXX, nCaXXI, nFeI, nFeII, 
           nFeIII, nFeIV, nFeV, nFeVI, nFeVII, nFeVIII, nFeIX, nFeX, nFeXI, nFeXII, 
           nFeXIII, nFeXIV, nFeXV, nFeXVI, nFeXVII, nFeXVIII, nFeXIX, nFeXX, nFeXXI, nFeXXII, 
           nFeXXIII, nFeXXIV, nFeXXV, nFeXXVI, nFeXXVII):

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
       + 1 * nCaII + 2 * nCaIII + 3 * nCaIV + 4 * nCaV + 5 * nCaVI
       + 6 * nCaVII + 7 * nCaVIII + 8 * nCaIX + 9 * nCaX + 10 * nCaXI + 11 * nCaXII
       + 12 * nCaXIII + 13 * nCaXIV + 14 * nCaXV + 15 * nCaXVI + 16 * nCaXVII + 17 * nCaXVIII
       + 18 * nCaXIX + 19 * nCaXX + 20 * nCaXXI + 1 * nFeII + 2 * nFeIII
       + 3 * nFeIV + 4 * nFeV + 5 * nFeVI + 6 * nFeVII + 7 * nFeVIII + 8 * nFeIX
       + 9 * nFeX + 10 * nFeXI + 11 * nFeXII + 12 * nFeXIII + 13 * nFeXIV + 14 * nFeXV
       + 15 * nFeXVI + 16 * nFeXVII + 17 * nFeXVIII + 18 * nFeXIX + 19 * nFeXX + 20 * nFeXXI
       + 21 * nFeXXII + 22 * nFeXXIII + 23 * nFeXXIV + 24 * nFeXXV + 25 * nFeXXVI + 26 * nFeXXVII
      
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
        + 10**cooling_rate_2d("FeII", T, ne, Temp_2d, elecDensity_2d) * nFeII * ne # cooling via FeII
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
        + 10**gCaI(Tx) * nCaI * ne
        + 10**gCaII(Tx) * nCaII * ne
        + 10**grain_recomb_cooling_rate * nCaII * ne # grain_recombination cooling!
        + 10**gCaIII(Tx) * nCaIII * ne
        + 10**grain_recomb_cooling_rate * nCaIII * ne # grain_recombination cooling!
        + 10**gCaIV(Tx) * nCaIV * ne
        + 10**gCaV(Tx) * nCaV * ne
        + 10**gCaVI(Tx) * nCaVI * ne
        + 10**gCaVII(Tx) * nCaVII * ne
        + 10**gCaVIII(Tx) * nCaVIII * ne
        + 10**gCaIX(Tx) * nCaIX * ne
        + 10**gCaX(Tx) * nCaX * ne
        + 10**gCaXI(Tx) * nCaXI * ne
        + 10**gCaXII(Tx) * nCaXII * ne
        + 10**gCaXIII(Tx) * nCaXIII * ne
        + 10**gCaXIV(Tx) * nCaXIV * ne
        + 10**gCaXV(Tx) * nCaXV * ne
        + 10**gCaXVI(Tx) * nCaXVI * ne
        + 10**gCaXVII(Tx) * nCaXVII * ne
        + 10**gCaXVIII(Tx) * nCaXVIII * ne
        + 10**gCaXIX(Tx) * nCaXIX * ne
        + 10**gCaXX(Tx) * nCaXX * ne
        + 10**gCaXXI(Tx) * nCaXXI * ne
        + 10**gFeI(Tx) * nFeI * ne
        + 10**grain_recomb_cooling_rate * nFeII * ne # grain_recombination cooling!
        + 10**gFeIII(Tx) * nFeIII * ne
        + 10**gFeIV(Tx) * nFeIV * ne
        + 10**gFeV(Tx) * nFeV * ne
        + 10**gFeVI(Tx) * nFeVI * ne
        + 10**gFeVII(Tx) * nFeVII * ne
        + 10**gFeVIII(Tx) * nFeVIII * ne
        + 10**gFeIX(Tx) * nFeIX * ne
        + 10**gFeX(Tx) * nFeX * ne
        + 10**gFeXI(Tx) * nFeXI * ne
        + 10**gFeXII(Tx) * nFeXII * ne
        + 10**gFeXIII(Tx) * nFeXIII * ne
        + 10**gFeXIV(Tx) * nFeXIV * ne
        + 10**gFeXV(Tx) * nFeXV * ne
        + 10**gFeXVI(Tx) * nFeXVI * ne
        + 10**gFeXVII(Tx) * nFeXVII * ne
        + 10**gFeXVIII(Tx) * nFeXVIII * ne
        + 10**gFeXIX(Tx) * nFeXIX * ne
        + 10**gFeXX(Tx) * nFeXX * ne
        + 10**gFeXXI(Tx) * nFeXXI * ne
        + 10**gFeXXII(Tx) * nFeXXII * ne
        + 10**gFeXXIII(Tx) * nFeXXIII * ne
        + 10**gFeXXIV(Tx) * nFeXXIV * ne
        + 10**gFeXXV(Tx) * nFeXXV * ne
        + 10**gFeXXVI(Tx) * nFeXXVI * ne
        + 10**gFeXXVII(Tx) * nFeXXVII * ne
        + LCompton)

  return Lamb



#----- func
def func(y):

  nHI, nHII, nHm, nHeI, nHeII, nHeIII, nCI, nCII, nCIII, nCIV, \
  nCV, nCVI, nCVII, nCm, nNI, nNII, nNIII, nNIV, nNV, nNVI, \
  nNVII, nNVIII, nOI, nOII, nOIII, nOIV, nOV, nOVI, nOVII, nOVIII, \
  nOIX, nOm, nNeI, nNeII, nNeIII, nNeIV, nNeV, nNeVI, nNeVII, nNeVIII, \
  nNeIX, nNeX, nNeXI, nMgI, nMgII, nMgIII, nMgIV, nMgV, nMgVI, nMgVII, \
  nMgVIII, nMgIX, nMgX, nMgXI, nMgXII, nMgXIII, nSiI, nSiII, nSiIII, nSiIV, \
  nSiV, nSiVI, nSiVII, nSiVIII, nSiIX, nSiX, nSiXI, nSiXII, nSiXIII, nSiXIV, \
  nSiXV, nSI, nSII, nSIII, nSIV, nSV, nSVI, nSVII, nSVIII, nSIX, \
  nSX, nSXI, nSXII, nSXIII, nSXIV, nSXV, nSXVI, nSXVII, nCaI, nCaII, \
  nCaIII, nCaIV, nCaV, nCaVI, nCaVII, nCaVIII, nCaIX, nCaX, nCaXI, nCaXII, \
  nCaXIII, nCaXIV, nCaXV, nCaXVI, nCaXVII, nCaXVIII, nCaXIX, nCaXX, nCaXXI, nFeI, \
  nFeII, nFeIII, nFeIV, nFeV, nFeVI, nFeVII, nFeVIII, nFeIX, nFeX, nFeXI, \
  nFeXII, nFeXIII, nFeXIV, nFeXV, nFeXVI, nFeXVII, nFeXVIII, nFeXIX, nFeXX, nFeXXI, \
  nFeXXII, nFeXXIII, nFeXXIV, nFeXXV, nFeXXVI, nFeXXVII, T = y

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
       + 16 * nSXVII + 1 * nCaII + 2 * nCaIII + 3 * nCaIV + 4 * nCaV + 5 * nCaVI
       + 6 * nCaVII + 7 * nCaVIII + 8 * nCaIX + 9 * nCaX + 10 * nCaXI + 11 * nCaXII
       + 12 * nCaXIII + 13 * nCaXIV + 14 * nCaXV + 15 * nCaXVI + 16 * nCaXVII + 17 * nCaXVIII
       + 18 * nCaXIX + 19 * nCaXX + 20 * nCaXXI + 1 * nFeII + 2 * nFeIII + 3 * nFeIV
       + 4 * nFeV + 5 * nFeVI + 6 * nFeVII + 7 * nFeVIII + 8 * nFeIX + 9 * nFeX
       + 10 * nFeXI + 11 * nFeXII + 12 * nFeXIII + 13 * nFeXIV + 14 * nFeXV + 15 * nFeXVI
       + 16 * nFeXVII + 17 * nFeXVIII + 18 * nFeXIX + 19 * nFeXX + 20 * nFeXXI + 21 * nFeXXII
       + 22 * nFeXXIII + 23 * nFeXXIV + 24 * nFeXXV + 25 * nFeXXVI + 26 * nFeXXVII
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
         + nSXIV + nSXV + nSXVI + nSXVII + nCaI + nCaII
         + nCaIII + nCaIV + nCaV + nCaVI + nCaVII + nCaVIII
         + nCaIX + nCaX + nCaXI + nCaXII + nCaXIII + nCaXIV
         + nCaXV + nCaXVI + nCaXVII + nCaXVIII + nCaXIX + nCaXX
         + nCaXXI + nFeI + nFeII + nFeIII + nFeIV + nFeV
         + nFeVI + nFeVII + nFeVIII + nFeIX + nFeX + nFeXI
         + nFeXII + nFeXIII + nFeXIV + nFeXV + nFeXVI + nFeXVII
         + nFeXVIII + nFeXIX + nFeXX + nFeXXI + nFeXXII + nFeXXIII
         + nFeXXIV + nFeXXV + nFeXXVI + nFeXXVII
         )

  grain_rec_HII_to_HI = grain_recomb_rate("HII", T, ne, G0, A_v, Temp, Psi)
  grain_rec_HeII_to_HeI = grain_recomb_rate("HeII", T, ne, G0, A_v, Temp, Psi)
  grain_rec_CII_to_CI = grain_recomb_rate("CII", T, ne, G0, A_v, Temp, Psi)
  grain_rec_OII_to_OI = grain_recomb_rate("OII", T, ne, G0, A_v, Temp, Psi)
  grain_rec_SiII_to_SiI = grain_recomb_rate("SiII", T, ne, G0, A_v, Temp, Psi)
  grain_rec_FeII_to_FeI = grain_recomb_rate("FeII", T, ne, G0, A_v, Temp, Psi)
  grain_rec_MgII_to_MgI = grain_recomb_rate("MgII", T, ne, G0, A_v, Temp, Psi)
  grain_rec_SII_to_SI = grain_recomb_rate("SII", T, ne, G0, A_v, Temp, Psi)
  grain_rec_CaII_to_CaI = grain_recomb_rate("CaII", T, ne, G0, A_v, Temp, Psi)
  grain_rec_CaIII_to_CaII = grain_recomb_rate("CaIII", T, ne, G0, A_v, Temp, Psi)

  const_CII_SiI_to_CI_SiII = 2.1000E-09 # constant/rates
  const_OI_e_to_Om_ = 1.5000E-15 # constant/rates
  const_CII_MgI_to_CI_MgII = 1.1000E-09 # constant/rates
  const_NII_MgI_to_NI_MgII = 1.2000E-09 # constant/rates
  const_CI_e_to_Cm_ = 2.2500E-15 # constant/rates
  const_SII_MgI_to_SI_MgII = 2.8000E-10 # constant/rates
  const_FeI_SiII_to_FeII_SiI = 1.9000E-10 # constant/rates
  const_FeI_CII_to_FeII_CI = 2.6000E-09 # constant/rates
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
             - 10**R_CaIV_to_CaIII_via_HI(Tx) * nCaIV * nHI
             - 10**R_CaV_to_CaIV_via_HI(Tx) * nCaV * nHI
             - 10**R_CaVI_to_CaV_via_HI(Tx) * nCaVI * nHI
             + 10**R_FeI_to_FeII_via_HII(Tx) * nFeI * nHII
             + 10**R_FeII_to_FeIII_via_HII(Tx) * nFeII * nHII
             - 10**R_FeIII_to_FeII_via_HI(Tx) * nFeIII * nHI
             - 10**R_FeIV_to_FeIII_via_HI(Tx) * nFeIV * nHI
             - 10**R_FeV_to_FeIV_via_HI(Tx) * nFeV * nHI
             - 10**R_FeVI_to_FeV_via_HI(Tx) * nFeVI * nHI
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
              + 10**R_CaIV_to_CaIII_via_HI(Tx) * nCaIV * nHI
              + 10**R_CaV_to_CaIV_via_HI(Tx) * nCaV * nHI
              + 10**R_CaVI_to_CaV_via_HI(Tx) * nCaVI * nHI
              - 10**R_FeI_to_FeII_via_HII(Tx) * nFeI * nHII
              - 10**R_FeII_to_FeIII_via_HII(Tx) * nFeII * nHII
              + 10**R_FeIII_to_FeII_via_HI(Tx) * nFeIII * nHI
              + 10**R_FeIV_to_FeIII_via_HI(Tx) * nFeIV * nHI
              + 10**R_FeV_to_FeIV_via_HI(Tx) * nFeV * nHI
              + 10**R_FeVI_to_FeV_via_HI(Tx) * nFeVI * nHI
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
             + const_FeI_CII_to_FeII_CI * nFeI * nCII # constant rate
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
              - const_FeI_CII_to_FeII_CI * nFeI * nCII # constant rate
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
              + const_FeI_SiII_to_FeII_SiI * nFeI * nSiII # constant rate
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
               - const_FeI_SiII_to_FeII_SiI * nFeI * nSiII # constant rate
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

  dnCaI_dt = (
              - 10**R_CaI_to_CaII_via_e(Tx) * nCaI * ne
              + 10**grain_rec_CaII_to_CaI * nCaII * ne # grain_recombination
              + 10**R_CaII_to_CaI_via_e(Tx) * nCaII * ne
             )

  dnCaII_dt = (
               + 10**R_CaI_to_CaII_via_e(Tx) * nCaI * ne
               - 10**R_CaII_to_CaI_via_e(Tx) * nCaII * ne
               - 10**grain_rec_CaII_to_CaI * nCaII * ne # grain_recombination
               + 10**grain_rec_CaIII_to_CaII * nCaIII * ne # grain_recombination
               - 10**R_CaII_to_CaIII_via_e(Tx) * nCaII * ne
               + 10**R_CaIII_to_CaII_via_e(Tx) * nCaIII * ne
              )

  dnCaIII_dt = (
                + 10**R_CaII_to_CaIII_via_e(Tx) * nCaII * ne
                - 10**R_CaIII_to_CaII_via_e(Tx) * nCaIII * ne
                - 10**grain_rec_CaIII_to_CaII * nCaIII * ne # grain_recombination
                - 10**R_CaIII_to_CaIV_via_e(Tx) * nCaIII * ne
                + 10**R_CaIV_to_CaIII_via_HI(Tx) * nCaIV * nHI
                + 10**R_CaIV_to_CaIII_via_e(Tx) * nCaIV * ne
               )

  dnCaIV_dt = (
               + 10**R_CaIII_to_CaIV_via_e(Tx) * nCaIII * ne
               - 10**R_CaIV_to_CaIII_via_HI(Tx) * nCaIV * nHI
               - 10**R_CaIV_to_CaIII_via_e(Tx) * nCaIV * ne
               - 10**R_CaIV_to_CaV_via_e(Tx) * nCaIV * ne
               + 10**R_CaV_to_CaIV_via_e(Tx) * nCaV * ne
               + 10**R_CaV_to_CaIV_via_HI(Tx) * nCaV * nHI
              )

  dnCaV_dt = (
              + 10**R_CaIV_to_CaV_via_e(Tx) * nCaIV * ne
              - 10**R_CaV_to_CaVI_via_e(Tx) * nCaV * ne
              - 10**R_CaV_to_CaIV_via_e(Tx) * nCaV * ne
              - 10**R_CaV_to_CaIV_via_HI(Tx) * nCaV * nHI
              + 10**R_CaVI_to_CaV_via_e(Tx) * nCaVI * ne
              + 10**R_CaVI_to_CaV_via_HI(Tx) * nCaVI * nHI
             )

  dnCaVI_dt = (
               + 10**R_CaV_to_CaVI_via_e(Tx) * nCaV * ne
               - 10**R_CaVI_to_CaVII_via_e(Tx) * nCaVI * ne
               - 10**R_CaVI_to_CaV_via_e(Tx) * nCaVI * ne
               - 10**R_CaVI_to_CaV_via_HI(Tx) * nCaVI * nHI
               + 10**R_CaVII_to_CaVI_via_e(Tx) * nCaVII * ne
              )

  dnCaVII_dt = (
                + 10**R_CaVI_to_CaVII_via_e(Tx) * nCaVI * ne
                - 10**R_CaVII_to_CaVI_via_e(Tx) * nCaVII * ne
                - 10**R_CaVII_to_CaVIII_via_e(Tx) * nCaVII * ne
                + 10**R_CaVIII_to_CaVII_via_e(Tx) * nCaVIII * ne
               )

  dnCaVIII_dt = (
                 + 10**R_CaVII_to_CaVIII_via_e(Tx) * nCaVII * ne
                 - 10**R_CaVIII_to_CaVII_via_e(Tx) * nCaVIII * ne
                 - 10**R_CaVIII_to_CaIX_via_e(Tx) * nCaVIII * ne
                 + 10**R_CaIX_to_CaVIII_via_e(Tx) * nCaIX * ne
                )

  dnCaIX_dt = (
               + 10**R_CaVIII_to_CaIX_via_e(Tx) * nCaVIII * ne
               - 10**R_CaIX_to_CaVIII_via_e(Tx) * nCaIX * ne
               - 10**R_CaIX_to_CaX_via_e(Tx) * nCaIX * ne
               + 10**R_CaX_to_CaIX_via_e(Tx) * nCaX * ne
              )

  dnCaX_dt = (
              + 10**R_CaIX_to_CaX_via_e(Tx) * nCaIX * ne
              - 10**R_CaX_to_CaIX_via_e(Tx) * nCaX * ne
              - 10**R_CaX_to_CaXI_via_e(Tx) * nCaX * ne
              + 10**R_CaXI_to_CaX_via_e(Tx) * nCaXI * ne
             )

  dnCaXI_dt = (
               + 10**R_CaX_to_CaXI_via_e(Tx) * nCaX * ne
               - 10**R_CaXI_to_CaX_via_e(Tx) * nCaXI * ne
               - 10**R_CaXI_to_CaXII_via_e(Tx) * nCaXI * ne
               + 10**R_CaXII_to_CaXI_via_e(Tx) * nCaXII * ne
              )

  dnCaXII_dt = (
                + 10**R_CaXI_to_CaXII_via_e(Tx) * nCaXI * ne
                - 10**R_CaXII_to_CaXI_via_e(Tx) * nCaXII * ne
                - 10**R_CaXII_to_CaXIII_via_e(Tx) * nCaXII * ne
                + 10**R_CaXIII_to_CaXII_via_e(Tx) * nCaXIII * ne
               )

  dnCaXIII_dt = (
                 + 10**R_CaXII_to_CaXIII_via_e(Tx) * nCaXII * ne
                 - 10**R_CaXIII_to_CaXII_via_e(Tx) * nCaXIII * ne
                 - 10**R_CaXIII_to_CaXIV_via_e(Tx) * nCaXIII * ne
                 + 10**R_CaXIV_to_CaXIII_via_e(Tx) * nCaXIV * ne
                )

  dnCaXIV_dt = (
                + 10**R_CaXIII_to_CaXIV_via_e(Tx) * nCaXIII * ne
                - 10**R_CaXIV_to_CaXIII_via_e(Tx) * nCaXIV * ne
                - 10**R_CaXIV_to_CaXV_via_e(Tx) * nCaXIV * ne
                + 10**R_CaXV_to_CaXIV_via_e(Tx) * nCaXV * ne
               )

  dnCaXV_dt = (
               + 10**R_CaXIV_to_CaXV_via_e(Tx) * nCaXIV * ne
               - 10**R_CaXV_to_CaXVI_via_e(Tx) * nCaXV * ne
               - 10**R_CaXV_to_CaXIV_via_e(Tx) * nCaXV * ne
               + 10**R_CaXVI_to_CaXV_via_e(Tx) * nCaXVI * ne
              )

  dnCaXVI_dt = (
                + 10**R_CaXV_to_CaXVI_via_e(Tx) * nCaXV * ne
                - 10**R_CaXVI_to_CaXV_via_e(Tx) * nCaXVI * ne
                - 10**R_CaXVI_to_CaXVII_via_e(Tx) * nCaXVI * ne
                + 10**R_CaXVII_to_CaXVI_via_e(Tx) * nCaXVII * ne
               )

  dnCaXVII_dt = (
                 + 10**R_CaXVI_to_CaXVII_via_e(Tx) * nCaXVI * ne
                 - 10**R_CaXVII_to_CaXVI_via_e(Tx) * nCaXVII * ne
                 - 10**R_CaXVII_to_CaXVIII_via_e(Tx) * nCaXVII * ne
                 + 10**R_CaXVIII_to_CaXVII_via_e(Tx) * nCaXVIII * ne
                )

  dnCaXVIII_dt = (
                  + 10**R_CaXVII_to_CaXVIII_via_e(Tx) * nCaXVII * ne
                  - 10**R_CaXVIII_to_CaXVII_via_e(Tx) * nCaXVIII * ne
                  - 10**R_CaXVIII_to_CaXIX_via_e(Tx) * nCaXVIII * ne
                  + 10**R_CaXIX_to_CaXVIII_via_e(Tx) * nCaXIX * ne
                 )

  dnCaXIX_dt = (
                + 10**R_CaXVIII_to_CaXIX_via_e(Tx) * nCaXVIII * ne
                - 10**R_CaXIX_to_CaXVIII_via_e(Tx) * nCaXIX * ne
                - 10**R_CaXIX_to_CaXX_via_e(Tx) * nCaXIX * ne
                + 10**R_CaXX_to_CaXIX_via_e(Tx) * nCaXX * ne
               )

  dnCaXX_dt = (
               + 10**R_CaXIX_to_CaXX_via_e(Tx) * nCaXIX * ne
               - 10**R_CaXX_to_CaXXI_via_e(Tx) * nCaXX * ne
               - 10**R_CaXX_to_CaXIX_via_e(Tx) * nCaXX * ne
               + 10**R_CaXXI_to_CaXX_via_e(Tx) * nCaXXI * ne
              )

  dnCaXXI_dt = (
                + 10**R_CaXX_to_CaXXI_via_e(Tx) * nCaXX * ne
                - 10**R_CaXXI_to_CaXX_via_e(Tx) * nCaXXI * ne
               )

  dnFeI_dt = (
              - 10**R_FeI_to_FeII_via_HII(Tx) * nFeI * nHII
              + 10**grain_rec_FeII_to_FeI * nFeII * ne # grain_recombination
              - 10**R_FeI_to_FeII_via_e(Tx) * nFeI * ne
              + 10**R_FeII_to_FeI_via_e(Tx) * nFeII * ne
              - const_FeI_SiII_to_FeII_SiI * nFeI * nSiII # constant rate
              - const_FeI_CII_to_FeII_CI * nFeI * nCII # constant rate
             )

  dnFeII_dt = (
               + 10**R_FeI_to_FeII_via_HII(Tx) * nFeI * nHII
               - 10**grain_rec_FeII_to_FeI * nFeII * ne # grain_recombination
               + 10**R_FeI_to_FeII_via_e(Tx) * nFeI * ne
               - 10**R_FeII_to_FeIII_via_HII(Tx) * nFeII * nHII
               - 10**R_FeII_to_FeI_via_e(Tx) * nFeII * ne
               - 10**R_FeII_to_FeIII_via_e(Tx) * nFeII * ne
               + 10**R_FeIII_to_FeII_via_e(Tx) * nFeIII * ne
               + 10**R_FeIII_to_FeII_via_HI(Tx) * nFeIII * nHI
               + const_FeI_SiII_to_FeII_SiI * nFeI * nSiII # constant rate
               + const_FeI_CII_to_FeII_CI * nFeI * nCII # constant rate
              )

  dnFeIII_dt = (
                + 10**R_FeII_to_FeIII_via_HII(Tx) * nFeII * nHII
                + 10**R_FeII_to_FeIII_via_e(Tx) * nFeII * ne
                - 10**R_FeIII_to_FeII_via_e(Tx) * nFeIII * ne
                - 10**R_FeIII_to_FeIV_via_e(Tx) * nFeIII * ne
                - 10**R_FeIII_to_FeII_via_HI(Tx) * nFeIII * nHI
                + 10**R_FeIV_to_FeIII_via_e(Tx) * nFeIV * ne
                + 10**R_FeIV_to_FeIII_via_HI(Tx) * nFeIV * nHI
               )

  dnFeIV_dt = (
               + 10**R_FeIII_to_FeIV_via_e(Tx) * nFeIII * ne
               - 10**R_FeIV_to_FeV_via_e(Tx) * nFeIV * ne
               - 10**R_FeIV_to_FeIII_via_e(Tx) * nFeIV * ne
               - 10**R_FeIV_to_FeIII_via_HI(Tx) * nFeIV * nHI
               + 10**R_FeV_to_FeIV_via_HI(Tx) * nFeV * nHI
               + 10**R_FeV_to_FeIV_via_e(Tx) * nFeV * ne
              )

  dnFeV_dt = (
              + 10**R_FeIV_to_FeV_via_e(Tx) * nFeIV * ne
              - 10**R_FeV_to_FeIV_via_HI(Tx) * nFeV * nHI
              - 10**R_FeV_to_FeIV_via_e(Tx) * nFeV * ne
              - 10**R_FeV_to_FeVI_via_e(Tx) * nFeV * ne
              + 10**R_FeVI_to_FeV_via_HI(Tx) * nFeVI * nHI
              + 10**R_FeVI_to_FeV_via_e(Tx) * nFeVI * ne
             )

  dnFeVI_dt = (
               + 10**R_FeV_to_FeVI_via_e(Tx) * nFeV * ne
               - 10**R_FeVI_to_FeV_via_HI(Tx) * nFeVI * nHI
               - 10**R_FeVI_to_FeV_via_e(Tx) * nFeVI * ne
               - 10**R_FeVI_to_FeVII_via_e(Tx) * nFeVI * ne
               + 10**R_FeVII_to_FeVI_via_e(Tx) * nFeVII * ne
              )

  dnFeVII_dt = (
                + 10**R_FeVI_to_FeVII_via_e(Tx) * nFeVI * ne
                - 10**R_FeVII_to_FeVI_via_e(Tx) * nFeVII * ne
                - 10**R_FeVII_to_FeVIII_via_e(Tx) * nFeVII * ne
                + 10**R_FeVIII_to_FeVII_via_e(Tx) * nFeVIII * ne
               )

  dnFeVIII_dt = (
                 + 10**R_FeVII_to_FeVIII_via_e(Tx) * nFeVII * ne
                 - 10**R_FeVIII_to_FeVII_via_e(Tx) * nFeVIII * ne
                 - 10**R_FeVIII_to_FeIX_via_e(Tx) * nFeVIII * ne
                 + 10**R_FeIX_to_FeVIII_via_e(Tx) * nFeIX * ne
                )

  dnFeIX_dt = (
               + 10**R_FeVIII_to_FeIX_via_e(Tx) * nFeVIII * ne
               - 10**R_FeIX_to_FeVIII_via_e(Tx) * nFeIX * ne
               - 10**R_FeIX_to_FeX_via_e(Tx) * nFeIX * ne
               + 10**R_FeX_to_FeIX_via_e(Tx) * nFeX * ne
              )

  dnFeX_dt = (
              + 10**R_FeIX_to_FeX_via_e(Tx) * nFeIX * ne
              - 10**R_FeX_to_FeIX_via_e(Tx) * nFeX * ne
              - 10**R_FeX_to_FeXI_via_e(Tx) * nFeX * ne
              + 10**R_FeXI_to_FeX_via_e(Tx) * nFeXI * ne
             )

  dnFeXI_dt = (
               + 10**R_FeX_to_FeXI_via_e(Tx) * nFeX * ne
               - 10**R_FeXI_to_FeX_via_e(Tx) * nFeXI * ne
               - 10**R_FeXI_to_FeXII_via_e(Tx) * nFeXI * ne
               + 10**R_FeXII_to_FeXI_via_e(Tx) * nFeXII * ne
              )

  dnFeXII_dt = (
                + 10**R_FeXI_to_FeXII_via_e(Tx) * nFeXI * ne
                - 10**R_FeXII_to_FeXI_via_e(Tx) * nFeXII * ne
                - 10**R_FeXII_to_FeXIII_via_e(Tx) * nFeXII * ne
                + 10**R_FeXIII_to_FeXII_via_e(Tx) * nFeXIII * ne
               )

  dnFeXIII_dt = (
                 + 10**R_FeXII_to_FeXIII_via_e(Tx) * nFeXII * ne
                 - 10**R_FeXIII_to_FeXII_via_e(Tx) * nFeXIII * ne
                 - 10**R_FeXIII_to_FeXIV_via_e(Tx) * nFeXIII * ne
                 + 10**R_FeXIV_to_FeXIII_via_e(Tx) * nFeXIV * ne
                )

  dnFeXIV_dt = (
                + 10**R_FeXIII_to_FeXIV_via_e(Tx) * nFeXIII * ne
                - 10**R_FeXIV_to_FeXIII_via_e(Tx) * nFeXIV * ne
                - 10**R_FeXIV_to_FeXV_via_e(Tx) * nFeXIV * ne
                + 10**R_FeXV_to_FeXIV_via_e(Tx) * nFeXV * ne
               )

  dnFeXV_dt = (
               + 10**R_FeXIV_to_FeXV_via_e(Tx) * nFeXIV * ne
               - 10**R_FeXV_to_FeXIV_via_e(Tx) * nFeXV * ne
               - 10**R_FeXV_to_FeXVI_via_e(Tx) * nFeXV * ne
               + 10**R_FeXVI_to_FeXV_via_e(Tx) * nFeXVI * ne
              )

  dnFeXVI_dt = (
                + 10**R_FeXV_to_FeXVI_via_e(Tx) * nFeXV * ne
                - 10**R_FeXVI_to_FeXV_via_e(Tx) * nFeXVI * ne
                - 10**R_FeXVI_to_FeXVII_via_e(Tx) * nFeXVI * ne
                + 10**R_FeXVII_to_FeXVI_via_e(Tx) * nFeXVII * ne
               )

  dnFeXVII_dt = (
                 + 10**R_FeXVI_to_FeXVII_via_e(Tx) * nFeXVI * ne
                 - 10**R_FeXVII_to_FeXVI_via_e(Tx) * nFeXVII * ne
                 - 10**R_FeXVII_to_FeXVIII_via_e(Tx) * nFeXVII * ne
                 + 10**R_FeXVIII_to_FeXVII_via_e(Tx) * nFeXVIII * ne
                )

  dnFeXVIII_dt = (
                  + 10**R_FeXVII_to_FeXVIII_via_e(Tx) * nFeXVII * ne
                  - 10**R_FeXVIII_to_FeXVII_via_e(Tx) * nFeXVIII * ne
                  - 10**R_FeXVIII_to_FeXIX_via_e(Tx) * nFeXVIII * ne
                  + 10**R_FeXIX_to_FeXVIII_via_e(Tx) * nFeXIX * ne
                 )

  dnFeXIX_dt = (
                + 10**R_FeXVIII_to_FeXIX_via_e(Tx) * nFeXVIII * ne
                - 10**R_FeXIX_to_FeXVIII_via_e(Tx) * nFeXIX * ne
                - 10**R_FeXIX_to_FeXX_via_e(Tx) * nFeXIX * ne
                + 10**R_FeXX_to_FeXIX_via_e(Tx) * nFeXX * ne
               )

  dnFeXX_dt = (
               + 10**R_FeXIX_to_FeXX_via_e(Tx) * nFeXIX * ne
               - 10**R_FeXX_to_FeXIX_via_e(Tx) * nFeXX * ne
               - 10**R_FeXX_to_FeXXI_via_e(Tx) * nFeXX * ne
               + 10**R_FeXXI_to_FeXX_via_e(Tx) * nFeXXI * ne
              )

  dnFeXXI_dt = (
                + 10**R_FeXX_to_FeXXI_via_e(Tx) * nFeXX * ne
                - 10**R_FeXXI_to_FeXX_via_e(Tx) * nFeXXI * ne
                - 10**R_FeXXI_to_FeXXII_via_e(Tx) * nFeXXI * ne
                + 10**R_FeXXII_to_FeXXI_via_e(Tx) * nFeXXII * ne
               )

  dnFeXXII_dt = (
                 + 10**R_FeXXI_to_FeXXII_via_e(Tx) * nFeXXI * ne
                 - 10**R_FeXXII_to_FeXXIII_via_e(Tx) * nFeXXII * ne
                 - 10**R_FeXXII_to_FeXXI_via_e(Tx) * nFeXXII * ne
                 + 10**R_FeXXIII_to_FeXXII_via_e(Tx) * nFeXXIII * ne
                )

  dnFeXXIII_dt = (
                  + 10**R_FeXXII_to_FeXXIII_via_e(Tx) * nFeXXII * ne
                  - 10**R_FeXXIII_to_FeXXIV_via_e(Tx) * nFeXXIII * ne
                  - 10**R_FeXXIII_to_FeXXII_via_e(Tx) * nFeXXIII * ne
                  + 10**R_FeXXIV_to_FeXXIII_via_e(Tx) * nFeXXIV * ne
                 )

  dnFeXXIV_dt = (
                 + 10**R_FeXXIII_to_FeXXIV_via_e(Tx) * nFeXXIII * ne
                 - 10**R_FeXXIV_to_FeXXV_via_e(Tx) * nFeXXIV * ne
                 - 10**R_FeXXIV_to_FeXXIII_via_e(Tx) * nFeXXIV * ne
                 + 10**R_FeXXV_to_FeXXIV_via_e(Tx) * nFeXXV * ne
                )

  dnFeXXV_dt = (
                + 10**R_FeXXIV_to_FeXXV_via_e(Tx) * nFeXXIV * ne
                - 10**R_FeXXV_to_FeXXIV_via_e(Tx) * nFeXXV * ne
                - 10**R_FeXXV_to_FeXXVI_via_e(Tx) * nFeXXV * ne
                + 10**R_FeXXVI_to_FeXXV_via_e(Tx) * nFeXXVI * ne
               )

  dnFeXXVI_dt = (
                 + 10**R_FeXXV_to_FeXXVI_via_e(Tx) * nFeXXV * ne
                 - 10**R_FeXXVI_to_FeXXV_via_e(Tx) * nFeXXVI * ne
                 - 10**R_FeXXVI_to_FeXXVII_via_e(Tx) * nFeXXVI * ne
                 + 10**R_FeXXVII_to_FeXXVI_via_e(Tx) * nFeXXVII * ne
                )

  dnFeXXVII_dt = (
                  + 10**R_FeXXVI_to_FeXXVII_via_e(Tx) * nFeXXVI * ne
                  - 10**R_FeXXVII_to_FeXXVI_via_e(Tx) * nFeXXVII * ne
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
                nSX, nSXI, nSXII, nSXIII, nSXIV, nSXV, nSXVI, nSXVII, nCaI, nCaII, 
                nCaIII, nCaIV, nCaV, nCaVI, nCaVII, nCaVIII, nCaIX, nCaX, nCaXI, nCaXII, 
                nCaXIII, nCaXIV, nCaXV, nCaXVI, nCaXVII, nCaXVIII, nCaXIX, nCaXX, nCaXXI, nFeI, 
                nFeII, nFeIII, nFeIV, nFeV, nFeVI, nFeVII, nFeVIII, nFeIX, nFeX, nFeXI, 
                nFeXII, nFeXIII, nFeXIV, nFeXV, nFeXVI, nFeXVII, nFeXVIII, nFeXIX, nFeXX, nFeXXI, 
                nFeXXII, nFeXXIII, nFeXXIV, nFeXXV, nFeXXVI, nFeXXVII)

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
           + 16 * dnSXVII_dt + 1 * dnCaII_dt + 2 * dnCaIII_dt + 3 * dnCaIV_dt + 4 * dnCaV_dt + 5 * dnCaVI_dt
           + 6 * dnCaVII_dt + 7 * dnCaVIII_dt + 8 * dnCaIX_dt + 9 * dnCaX_dt + 10 * dnCaXI_dt + 11 * dnCaXII_dt
           + 12 * dnCaXIII_dt + 13 * dnCaXIV_dt + 14 * dnCaXV_dt + 15 * dnCaXVI_dt + 16 * dnCaXVII_dt + 17 * dnCaXVIII_dt
           + 18 * dnCaXIX_dt + 19 * dnCaXX_dt + 20 * dnCaXXI_dt + 1 * dnFeII_dt + 2 * dnFeIII_dt + 3 * dnFeIV_dt
           + 4 * dnFeV_dt + 5 * dnFeVI_dt + 6 * dnFeVII_dt + 7 * dnFeVIII_dt + 8 * dnFeIX_dt + 9 * dnFeX_dt
           + 10 * dnFeXI_dt + 11 * dnFeXII_dt + 12 * dnFeXIII_dt + 13 * dnFeXIV_dt + 14 * dnFeXV_dt + 15 * dnFeXVI_dt
           + 16 * dnFeXVII_dt + 17 * dnFeXVIII_dt + 18 * dnFeXIX_dt + 19 * dnFeXX_dt + 20 * dnFeXXI_dt + 21 * dnFeXXII_dt
           + 22 * dnFeXXIII_dt + 23 * dnFeXXIV_dt + 24 * dnFeXXV_dt + 25 * dnFeXXVI_dt + 26 * dnFeXXVII_dt
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
              + dnSXIV_dt + dnSXV_dt + dnSXVI_dt + dnSXVII_dt + dnCaI_dt + dnCaII_dt
              + dnCaIII_dt + dnCaIV_dt + dnCaV_dt + dnCaVI_dt + dnCaVII_dt + dnCaVIII_dt
              + dnCaIX_dt + dnCaX_dt + dnCaXI_dt + dnCaXII_dt + dnCaXIII_dt + dnCaXIV_dt
              + dnCaXV_dt + dnCaXVI_dt + dnCaXVII_dt + dnCaXVIII_dt + dnCaXIX_dt + dnCaXX_dt
              + dnCaXXI_dt + dnFeI_dt + dnFeII_dt + dnFeIII_dt + dnFeIV_dt + dnFeV_dt
              + dnFeVI_dt + dnFeVII_dt + dnFeVIII_dt + dnFeIX_dt + dnFeX_dt + dnFeXI_dt
              + dnFeXII_dt + dnFeXIII_dt + dnFeXIV_dt + dnFeXV_dt + dnFeXVI_dt + dnFeXVII_dt
              + dnFeXVIII_dt + dnFeXIX_dt + dnFeXX_dt + dnFeXXI_dt + dnFeXXII_dt + dnFeXXIII_dt
              + dnFeXXIV_dt + dnFeXXV_dt + dnFeXXVI_dt + dnFeXXVII_dt
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
          dnSXIV_dt, dnSXV_dt, dnSXVI_dt, dnSXVII_dt, dnCaI_dt, dnCaII_dt,
          dnCaIII_dt, dnCaIV_dt, dnCaV_dt, dnCaVI_dt, dnCaVII_dt, dnCaVIII_dt,
          dnCaIX_dt, dnCaX_dt, dnCaXI_dt, dnCaXII_dt, dnCaXIII_dt, dnCaXIV_dt,
          dnCaXV_dt, dnCaXVI_dt, dnCaXVII_dt, dnCaXVIII_dt, dnCaXIX_dt, dnCaXX_dt,
          dnCaXXI_dt, dnFeI_dt, dnFeII_dt, dnFeIII_dt, dnFeIV_dt, dnFeV_dt,
          dnFeVI_dt, dnFeVII_dt, dnFeVIII_dt, dnFeIX_dt, dnFeX_dt, dnFeXI_dt,
          dnFeXII_dt, dnFeXIII_dt, dnFeXIV_dt, dnFeXV_dt, dnFeXVI_dt, dnFeXVII_dt,
          dnFeXVIII_dt, dnFeXIX_dt, dnFeXX_dt, dnFeXXI_dt, dnFeXXII_dt, dnFeXXIII_dt,
          dnFeXXIV_dt, dnFeXXV_dt, dnFeXXVI_dt, dnFeXXVII_dt, dT_dt
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

Ca_solar = 10**(-5.64)
nCa = Ca_solar * nH

Fe_solar = 10**(-4.55)
nFe = Fe_solar * nH

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

nCa0_i = 1e-20 * nCa
nCa1_i = 1e-19 * nCa
nCa2_i = 1e-18 * nCa
nCa3_i = 1e-17 * nCa
nCa4_i = 1e-16 * nCa
nCa5_i = 1e-15 * nCa
nCa6_i = 1e-14 * nCa
nCa7_i = 1e-13 * nCa
nCa8_i = 1e-12 * nCa
nCa9_i = 1e-11 * nCa
nCa10_i = 1e-10 * nCa
nCa11_i = 1e-9 * nCa
nCa12_i = 1e-8 * nCa
nCa13_i = 1e-7 * nCa
nCa14_i = 1e-6 * nCa
nCa15_i = 1e-5 * nCa
nCa16_i = 1e-4 * nCa
nCa17_i = 1e-3 * nCa
nCa18_i = 1e-2 * nCa
nCa19_i = 1e-2 * nCa
nCa20_i = nCa - (nCa0_i + nCa1_i + nCa2_i + nCa3_i + nCa4_i + nCa5_i + nCa6_i + nCa7_i + nCa8_i + nCa9_i + nCa10_i + nCa11_i + nCa12_i + nCa13_i + nCa14_i + nCa15_i + nCa16_i + nCa17_i + nCa18_i + nCa19_i)

nFe0_i = 1e-26 * nFe
nFe1_i = 1e-25 * nFe
nFe2_i = 1e-24 * nFe
nFe3_i = 1e-23 * nFe
nFe4_i = 1e-22 * nFe
nFe5_i = 1e-21 * nFe
nFe6_i = 1e-20 * nFe
nFe7_i = 1e-19 * nFe
nFe8_i = 1e-18 * nFe
nFe9_i = 1e-17 * nFe
nFe10_i = 1e-16 * nFe
nFe11_i = 1e-15 * nFe
nFe12_i = 1e-14 * nFe
nFe13_i = 1e-13 * nFe
nFe14_i = 1e-12 * nFe
nFe15_i = 1e-11 * nFe
nFe16_i = 1e-10 * nFe
nFe17_i = 1e-9 * nFe
nFe18_i = 1e-8 * nFe
nFe19_i = 1e-7 * nFe
nFe20_i = 1e-6 * nFe
nFe21_i = 1e-5 * nFe
nFe22_i = 1e-4 * nFe
nFe23_i = 1e-3 * nFe
nFe24_i = 1e-2 * nFe
nFe25_i = 1e-2 * nFe
nFe26_i = nFe - (nFe0_i + nFe1_i + nFe2_i + nFe3_i + nFe4_i + nFe5_i + nFe6_i + nFe7_i + nFe8_i + nFe9_i + nFe10_i + nFe11_i + nFe12_i + nFe13_i + nFe14_i + nFe15_i + nFe16_i + nFe17_i + nFe18_i + nFe19_i + nFe20_i + nFe21_i + nFe22_i + nFe23_i + nFe24_i + nFe25_i)


T_i = 10**6.00 # initial Temperature

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
      nS13_i, nS14_i, nS15_i, nS16_i, nCa0_i, nCa1_i, 
      nCa2_i, nCa3_i, nCa4_i, nCa5_i, nCa6_i, nCa7_i, 
      nCa8_i, nCa9_i, nCa10_i, nCa11_i, nCa12_i, nCa13_i, 
      nCa14_i, nCa15_i, nCa16_i, nCa17_i, nCa18_i, nCa19_i, 
      nCa20_i, nFe0_i, nFe1_i, nFe2_i, nFe3_i, nFe4_i, 
      nFe5_i, nFe6_i, nFe7_i, nFe8_i, nFe9_i, nFe10_i, 
      nFe11_i, nFe12_i, nFe13_i, nFe14_i, nFe15_i, nFe16_i, 
      nFe17_i, nFe18_i, nFe19_i, nFe20_i, nFe21_i, nFe22_i, 
      nFe23_i, nFe24_i, nFe25_i, nFe26_i, 
      T_i
      ]

A_v = 1.0
G0 = 0.01
dust_ratio = 0.01

#t_span = np.arange(1*3.16e7, 100*3.16e7, 3.16e7)
#t = t_span[0]

TA = time.time()
dy = func(y0)
print('Elapsed time = ', time.time() - TA)

dt = 1000

res = []

k = 0

t = 0

while k < 10000:
  for i, (y0x, dyx) in enumerate(zip(y0, dy)):
    y0[i] = y0x + dyx * dt
    y0 = np.abs(y0)
    #print('np.sort(y0) = ', np.sort(y0))
  dy = func(y0)
  
  res.append([t, y0[-1]])
  
  k += 1
  
  t += dt

res = np.array(res)
t = res[:, 0]
T = res[:, 1]

print()
print('np.log10(T) = ', np.log10(T))

t_yrs = t / 3.16e7

plt.scatter(t_yrs, np.log10(T), s = 1)
plt.ylim(4, 7)

plt.show()










