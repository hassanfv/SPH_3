#----- Lambda
def Lambda(T, nHI, nHII, nHm, nHeI, nHeII, nHeIII, nOI, nOII, nOIII, nOIV, nOV, 
           nOVI, nOVII, nOVIII, nOIX, nOm):

  Tx = np.log10(T)

  ne = (
         1 * nHII - nHm + (nHeII + 2.0 * nHeIII) + 1 * nOII + 2 * nOIII + 3 * nOIV
       + 4 * nOV + 5 * nOVI + 6 * nOVII + 7 * nOVIII + 8 * nOIX - 1 * nOm
       )

  cFree = (
            1 * nHII + nHeII + 4.0 * nHeIII + 1 * nOII + 4 * nOIII + 9 * nOIV
          + 16 * nOV + 25 * nOVI + 36 * nOVII + 49 * nOVIII + 64 * nOIX )

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
        + 10**cooling_rate_4d("OI", T, nHI, ne, nHII, Temp_4d, HIDensity_4d, elecDensity_4d, HIIDensity_4d) * nOI * ne # cooling via OI
        + 10**gOII(Tx) * nOII * ne
        + 10**grain_recomb_cooling_rate * nOII * ne # grain_recombination cooling!
        + 10**gOIII(Tx) * nOIII * ne
        + 10**gOIV(Tx) * nOIV * ne
        + 10**gOV(Tx) * nOV * ne
        + 10**gOVI(Tx) * nOVI * ne
        + 10**gOVII(Tx) * nOVII * ne
        + 10**gOVIII(Tx) * nOVIII * ne
        + 10**gOIX(Tx) * nOIX * ne
        + gfree(T) * ne * cFree # free-free emission
        + LCompton)

  return Lamb