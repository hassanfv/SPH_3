#----- func
def func(t, y):

  nHI, nHII, nHm, nHeI, nHeII, nHeIII, nNI, nNII, nNIII, nNIV, \
  nNV, nNVI, nNVII, nNVIII, T = y

  Tx = np.log10(T)

  ne = (
       + 1 * nHII + -1 * nHm + 1 * nHeII + 2 * nHeIII + 1 * nNII + 2 * nNIII
       + 3 * nNIV + 4 * nNV + 5 * nNVI + 6 * nNVII + 7 * nNVIII
       )

  ntot = (
           ne + nHI + nHII + nHm + nHeI + nHeII + nHeIII
         + nNI + nNII + nNIII + nNIV + nNV + nNVI
         + nNVII + nNVIII
         )

  grain_rec_HII_to_HI = grain_recomb_rate("HII", T, ne, G0, A_v, Temp, Psi)
  grain_rec_HeII_to_HeI = grain_recomb_rate("HeII", T, ne, G0, A_v, Temp, Psi)


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
                T, nHI, nHII, nHm, nHeI, nHeII, nHeIII, nNI, nNII, nNIII, nNIV, 
                nNV, nNVI, nNVII, nNVIII)

  dne_dt = (
           + 1 * dnHII_dt + -1 * dnHm_dt + 1 * dnHeII_dt + 2 * dnHeIII_dt + 1 * dnNII_dt + 2 * dnNIII_dt
           + 3 * dnNIV_dt + 4 * dnNV_dt + 5 * dnNVI_dt + 6 * dnNVII_dt + 7 * dnNVIII_dt
           )

  dntot_dt = (
                dne_dt + dnHI_dt + dnHII_dt + dnHm_dt + dnHeI_dt + dnHeII_dt + dnHeIII_dt
              + dnNI_dt + dnNII_dt + dnNIII_dt + dnNIV_dt + dnNV_dt + dnNVI_dt
              + dnNVII_dt + dnNVIII_dt
             )

  dT_dt = -1.0 * (gamma - 1.0) / kB / ntot * (Lamb + 1. / (gamma - 1.) * kB * T * dntot_dt)

  return [
          dnHI_dt, dnHII_dt, dnHm_dt, dnHeI_dt, dnHeII_dt, dnHeIII_dt,
          dnNI_dt, dnNII_dt, dnNIII_dt, dnNIV_dt, dnNV_dt, dnNVI_dt,
          dnNVII_dt, dnNVIII_dt, dT_dt
         ]

