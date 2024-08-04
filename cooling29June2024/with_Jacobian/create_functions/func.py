#----- func
def func(t, y):

  nHI, nHII, nHm, nHeI, nHeII, nHeIII, nCI, nCII, nCIII, nCIV, \
  nCV, nCVI, nCVII, nCm, nNI, nNII, nNIII, nNIV, nNV, nNVI, \
  nNVII, nNVIII, nOI, nOII, nOIII, nOIV, nOV, nOVI, nOVII, nOVIII, \
  nOIX, nOm, T = y

  Tx = np.log10(T)

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

