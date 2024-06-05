#----- Lambda
def Lambda(T, nHI, nHII, nHeI, nHeII, nHeIII, nCI, nCII, nCIII, nCIV, nCV, 
           nCVI, nCVII):

  Tx = np.log10(T)

  ne = (nHII + (nHeII + 2.0 * nHeIII) + 1 * nCII + 2 * nCIII + 3 * nCIV
     + 4 * nCV + 5 * nCVI + 6 * nCVII)

  cFree = (nHII + nHeII + 4.0 * nHeIII + 1 * nCII + 4 * nCIII + 9 * nCIV
        + 16 * nCV + 25 * nCVI + 36 * nCVII)

  #----- # Glover & Jappsen - 2007 -----
  z = 0.0 # current time redshift!
  TCMB_0 = 2.7255
  TCMB = TCMB_0 * (1.0 + z)
  LCompton = 1.017e-37 * TCMB**4 * (T - TCMB) * ne
  #--------------------------------------

  grain_cool_H_1_0 = grain_cool_rate(T, ne, nHI, nHII, G0, A_v, dust_ratio, Temp, Psi)
  grain_cool_He_1_0 = grain_cool_rate(T, ne, nHI, nHII, G0, A_v, dust_ratio, Temp, Psi)
  grain_cool_C_1_0 = grain_cool_rate(T, ne, nHI, nHII, G0, A_v, dust_ratio, Temp, Psi)

  Lamb = (
          10**g1(Tx) * ne * nHI  # H0
        + 10**g2(Tx) * ne * nHII # Hp
        + 10**grain_cool_H_1_0 * nHII * ne
        + 10**g3(Tx) * nHeI * ne # He0
        + 10**g4(Tx) * nHeII * ne # Hep
        + 10**grain_cool_He_1_0 * nHeII * ne
        + 10**g5(Tx) * nHeIII * ne# Hep
        + 10**CI_cooling_rate(T, nHI, ne, nHII, Temp_4d, HIDensity_4d, elecDensity_4d, HIIDensity_4d) * nCI * ne # cooling via CI
        + 10**CII_cooling_rate(T, ne, Temp_2d, elecDensity_2d) * nCII * ne # cooling via CII
        + 10**grain_cool_C_1_0 * nCII * ne
        + 10**gCIII(Tx) * nCIII * ne
        + 10**gCIV(Tx) * nCIV * ne
        + 10**gCV(Tx) * nCV * ne
        + 10**gCVI(Tx) * nCVI * ne
        + 10**gCVII(Tx) * nCVII * ne
        + gfree(T) * ne * cFree # free-free emission
        + LCompton)

  return Lamb