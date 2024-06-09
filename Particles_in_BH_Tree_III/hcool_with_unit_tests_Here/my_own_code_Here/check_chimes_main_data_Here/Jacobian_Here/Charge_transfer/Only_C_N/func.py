def func(t, y):

  nHI, nHII, nHm, nHeI, nHeII, nHeIII, nCI, nCII, nCIII, nCIV, 
  nCV, nCVI, nCVII, nCm, T = y

  Tx = np.log10(T)

  ne = (
       + 1 * nHII + -1 * nHm + 1 * nHeII + 2 * nHeIII + 1 * nCII + 2 * nCIII
       + 3 * nCIV + 4 * nCV + 5 * nCVI + 6 * nCVII + -1 * nCm
       )

  ntot = (
           ne + nHI + nHII + nHm + nHeI + nHeII + nHeIII
         + nCI + nCII + nCIII + nCIV + nCV + nCVI
         + nCVII + nCm
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

  Lamb = Lambda(
                T, nHI, nHII, nHm, nHeI, nHeII, nHeIII, nCI, nCII, nCIII, nCIV, 
                nCV, nCVI, nCVII, nCm)

  dne_dt = (
           + 1 * dnHII_dt + -1 * dnHm_dt + 1 * dnHeII_dt + 2 * dnHeIII_dt + 1 * dnCII_dt + 2 * dnCIII_dt
           + 3 * dnCIV_dt + 4 * dnCV_dt + 5 * dnCVI_dt + 6 * dnCVII_dt + -1 * dnCm_dt
           )

