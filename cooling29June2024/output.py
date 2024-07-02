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