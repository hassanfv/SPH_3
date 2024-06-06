dnCI_dt = (
           - 10**R_CI_to_CII_via_HeII(Tx) * nCI * nHeII
           + 10**grain_rec_CII_to_CI * nCII * ne) # grain_recombination
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
            - 10**grain_rec_CII_to_CI * nCII * ne) # grain_recombination
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
           + 10**grain_rec_OII_to_OI * nOII * ne) # grain_recombination
           - 10**R_OI_to_OII_via_e(Tx) * nOI * ne
           - 10**R_OI_to_OII_via_HII(Tx) * nOI * nHII
           + 10**R_OII_to_OI_via_HI(Tx) * nOII * nHI
           + 10**R_OII_to_OI_via_e(Tx) * nOII * ne
           + 10**R_Om_to_OI_via_HII(Tx) * nOm * nHII
           - const_OI_e_to_Om_ * nOI * ne # constant rate
          )

dnOII_dt = (
            + 10**R_OI_to_OII_via_HeII(Tx) * nOI * nHeII
            - 10**grain_rec_OII_to_OI * nOII * ne) # grain_recombination
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
            + 10**grain_rec_MgII_to_MgI * nMgII * ne) # grain_recombination
            - 10**R_MgI_to_MgII_via_e(Tx) * nMgI * ne
            + 10**R_MgII_to_MgI_via_e(Tx) * nMgII * ne
            - const_CII_MgI_to_CI_MgII * nCII * nMgI # constant rate
            - const_NII_MgI_to_NI_MgII * nNII * nMgI # constant rate
            - const_SII_MgI_to_SI_MgII * nSII * nMgI # constant rate
            - const_SiII_MgI_to_SiI_MgII * nSiII * nMgI # constant rate
           )

dnMgII_dt = (
             + 10**R_MgI_to_MgII_via_HII(Tx) * nMgI * nHII
             - 10**grain_rec_MgII_to_MgI * nMgII * ne) # grain_recombination
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
            + 10**grain_rec_SiII_to_SiI * nSiII * ne) # grain_recombination
            - 10**R_SiI_to_SiII_via_HII(Tx) * nSiI * nHII
            - 10**R_SiI_to_SiII_via_e(Tx) * nSiI * ne
            + 10**R_SiII_to_SiI_via_e(Tx) * nSiII * ne
            - const_CII_SiI_to_CI_SiII * nCII * nSiI # constant rate
            + const_FeI_SiII_to_FeII_SiI * nFeI * nSiII # constant rate
            + const_SiII_MgI_to_SiI_MgII * nSiII * nMgI # constant rate
           )

dnSiII_dt = (
             + 10**R_SiI_to_SiII_via_HeII(Tx) * nSiI * nHeII
             - 10**grain_rec_SiII_to_SiI * nSiII * ne) # grain_recombination
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
           + 10**grain_rec_SII_to_SI * nSII * ne) # grain_recombination
           - 10**R_SI_to_SII_via_e(Tx) * nSI * ne
           + 10**R_SII_to_SI_via_HI(Tx) * nSII * nHI
           + 10**R_SII_to_SI_via_e(Tx) * nSII * ne
           + const_SII_MgI_to_SI_MgII * nSII * nMgI # constant rate
          )

dnSII_dt = (
            + 10**R_SI_to_SII_via_HII(Tx) * nSI * nHII
            - 10**grain_rec_SII_to_SI * nSII * ne) # grain_recombination
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
            + 10**grain_rec_CaII_to_CaI * nCaII * ne) # grain_recombination
            + 10**R_CaII_to_CaI_via_e(Tx) * nCaII * ne
           )

dnCaII_dt = (
             + 10**R_CaI_to_CaII_via_e(Tx) * nCaI * ne
             - 10**R_CaII_to_CaI_via_e(Tx) * nCaII * ne
             - 10**grain_rec_CaII_to_CaI * nCaII * ne) # grain_recombination
             + 10**grain_rec_CaIII_to_CaII * nCaIII * ne) # grain_recombination
             - 10**R_CaII_to_CaIII_via_e(Tx) * nCaII * ne
             + 10**R_CaIII_to_CaII_via_e(Tx) * nCaIII * ne
            )

dnCaIII_dt = (
              + 10**R_CaII_to_CaIII_via_e(Tx) * nCaII * ne
              - 10**R_CaIII_to_CaII_via_e(Tx) * nCaIII * ne
              - 10**grain_rec_CaIII_to_CaII * nCaIII * ne) # grain_recombination
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
            + 10**grain_rec_FeII_to_FeI * nFeII * ne) # grain_recombination
            - 10**R_FeI_to_FeII_via_e(Tx) * nFeI * ne
            + 10**R_FeII_to_FeI_via_e(Tx) * nFeII * ne
            - const_FeI_SiII_to_FeII_SiI * nFeI * nSiII # constant rate
            - const_FeI_CII_to_FeII_CI * nFeI * nCII # constant rate
           )

dnFeII_dt = (
             + 10**R_FeI_to_FeII_via_HII(Tx) * nFeI * nHII
             - 10**grain_rec_FeII_to_FeI * nFeII * ne) # grain_recombination
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

