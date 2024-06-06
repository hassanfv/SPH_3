dnCI_dt = (
           - 10**R_CI_to_CII_via_HeII(Tx) * nCI * nHeII
           + 10**grain_rec_CII_to_CI * nCII * ne) # grain_recombination
           - 10**R_CI_to_CII_via_HII(Tx) * nCI * nHII
           - 10**R_CI_to_CII_via_e(Tx) * nCI * ne
           + 10**R_CII_to_CI_via_HI(Tx) * nCII * nHI
           + 10**R_CII_to_CI_via_e(Tx) * nCII * ne
           + 10**R_Cm_to_CI_via_HII(Tx) * nCm * nHII
           - const_CI_e_to_Cm_ * nCI * ne # constant rate
          )

dnCII_dt = (
            + 10**R_CI_to_CII_via_HeII(Tx) * nCI * nHeII
            - 10**grain_rec_CII_to_CI * nCII * ne) # grain_recombination
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

