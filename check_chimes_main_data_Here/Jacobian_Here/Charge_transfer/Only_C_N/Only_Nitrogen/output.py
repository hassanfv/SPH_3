#----- Lambda
def Lambda(T, nHI, nHII, nHm, nHeI, nHeII, nHeIII, nNI, nNII, nNIII, nNIV, nNV, 
           nNVI, nNVII, nNVIII):

  Tx = np.log10(T)

  ne = (
         1 * nHII - nHm + (nHeII + 2.0 * nHeIII) + 1 * nNII + 2 * nNIII + 3 * nNIV
       + 4 * nNV + 5 * nNVI + 6 * nNVII + 7 * nNVIII
       )

  cFree = (
            1 * nHII + nHeII + 4.0 * nHeIII + 1 * nNII + 4 * nNIII + 9 * nNIV
          + 16 * nNV + 25 * nNVI + 36 * nNVII + 49 * nNVIII )

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
        + 10**cooling_rate_2d("NII", T, ne, Temp_2d, elecDensity_2d) * nNII * ne # cooling via NII
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