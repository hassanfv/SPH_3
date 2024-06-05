dnMgI_dt = (
            - 10**R_MgI_to_MgII_via_HII(Tx) * nMgI * nHII
            + 10**grain_rec_MgII_to_MgI * nMgII * ne) # grain_recombination
            - 10**R_MgI_to_MgII_via_e(Tx) * nMgI * ne
            + 10**R_MgII_to_MgI_via_e(Tx) * nMgII * ne)

dnMgII_dt = (
             + 10**R_MgI_to_MgII_via_HII(Tx) * nMgI * nHII
             - 10**grain_rec_MgII_to_MgI * nMgII * ne) # grain_recombination
             + 10**R_MgI_to_MgII_via_e(Tx) * nMgI * ne
             - 10**R_MgII_to_MgIII_via_HII(Tx) * nMgII * nHII
             - 10**R_MgII_to_MgI_via_e(Tx) * nMgII * ne
             - 10**R_MgII_to_MgIII_via_e(Tx) * nMgII * ne
             + 10**R_MgIII_to_MgII_via_HI(Tx) * nMgIII * nHI
             + 10**R_MgIII_to_MgII_via_e(Tx) * nMgIII * ne)

dnMgIII_dt = (
              + 10**R_MgII_to_MgIII_via_HII(Tx) * nMgII * nHII
              + 10**R_MgII_to_MgIII_via_e(Tx) * nMgII * ne
              - 10**R_MgIII_to_MgII_via_HI(Tx) * nMgIII * nHI
              - 10**R_MgIII_to_MgII_via_e(Tx) * nMgIII * ne
              - 10**R_MgIII_to_MgIV_via_e(Tx) * nMgIII * ne
              + 10**R_MgIV_to_MgIII_via_HeI(Tx) * nMgIV * nHeI
              + 10**R_MgIV_to_MgIII_via_e(Tx) * nMgIV * ne
              + 10**R_MgIV_to_MgIII_via_HI(Tx) * nMgIV * nHI)

dnMgIV_dt = (
             + 10**R_MgIII_to_MgIV_via_e(Tx) * nMgIII * ne
             - 10**R_MgIV_to_MgIII_via_HeI(Tx) * nMgIV * nHeI
             - 10**R_MgIV_to_MgIII_via_e(Tx) * nMgIV * ne
             - 10**R_MgIV_to_MgV_via_e(Tx) * nMgIV * ne
             - 10**R_MgIV_to_MgIII_via_HI(Tx) * nMgIV * nHI
             + 10**R_MgV_to_MgIV_via_HeI(Tx) * nMgV * nHeI
             + 10**R_MgV_to_MgIV_via_HI(Tx) * nMgV * nHI
             + 10**R_MgV_to_MgIV_via_e(Tx) * nMgV * ne)

dnMgV_dt = (
            + 10**R_MgIV_to_MgV_via_e(Tx) * nMgIV * ne
            - 10**R_MgV_to_MgIV_via_HeI(Tx) * nMgV * nHeI
            - 10**R_MgV_to_MgIV_via_HI(Tx) * nMgV * nHI
            - 10**R_MgV_to_MgIV_via_e(Tx) * nMgV * ne
            - 10**R_MgV_to_MgVI_via_e(Tx) * nMgV * ne
            + 10**R_MgVI_to_MgV_via_HI(Tx) * nMgVI * nHI
            + 10**R_MgVI_to_MgV_via_e(Tx) * nMgVI * ne)

dnMgVI_dt = (
             + 10**R_MgV_to_MgVI_via_e(Tx) * nMgV * ne
             - 10**R_MgVI_to_MgV_via_HI(Tx) * nMgVI * nHI
             - 10**R_MgVI_to_MgVII_via_e(Tx) * nMgVI * ne
             - 10**R_MgVI_to_MgV_via_e(Tx) * nMgVI * ne
             + 10**R_MgVII_to_MgVI_via_e(Tx) * nMgVII * ne)

dnMgVII_dt = (
              + 10**R_MgVI_to_MgVII_via_e(Tx) * nMgVI * ne
              - 10**R_MgVII_to_MgVIII_via_e(Tx) * nMgVII * ne
              - 10**R_MgVII_to_MgVI_via_e(Tx) * nMgVII * ne
              + 10**R_MgVIII_to_MgVII_via_e(Tx) * nMgVIII * ne)

dnMgVIII_dt = (
               + 10**R_MgVII_to_MgVIII_via_e(Tx) * nMgVII * ne
               - 10**R_MgVIII_to_MgVII_via_e(Tx) * nMgVIII * ne
               - 10**R_MgVIII_to_MgIX_via_e(Tx) * nMgVIII * ne
               + 10**R_MgIX_to_MgVIII_via_e(Tx) * nMgIX * ne)

dnMgIX_dt = (
             + 10**R_MgVIII_to_MgIX_via_e(Tx) * nMgVIII * ne
             - 10**R_MgIX_to_MgVIII_via_e(Tx) * nMgIX * ne
             - 10**R_MgIX_to_MgX_via_e(Tx) * nMgIX * ne
             + 10**R_MgX_to_MgIX_via_e(Tx) * nMgX * ne)

dnMgX_dt = (
            + 10**R_MgIX_to_MgX_via_e(Tx) * nMgIX * ne
            - 10**R_MgX_to_MgIX_via_e(Tx) * nMgX * ne
            - 10**R_MgX_to_MgXI_via_e(Tx) * nMgX * ne
            + 10**R_MgXI_to_MgX_via_e(Tx) * nMgXI * ne)

dnMgXI_dt = (
             + 10**R_MgX_to_MgXI_via_e(Tx) * nMgX * ne
             - 10**R_MgXI_to_MgX_via_e(Tx) * nMgXI * ne
             - 10**R_MgXI_to_MgXII_via_e(Tx) * nMgXI * ne
             + 10**R_MgXII_to_MgXI_via_e(Tx) * nMgXII * ne)

dnMgXII_dt = (
              + 10**R_MgXI_to_MgXII_via_e(Tx) * nMgXI * ne
              - 10**R_MgXII_to_MgXIII_via_e(Tx) * nMgXII * ne
              - 10**R_MgXII_to_MgXI_via_e(Tx) * nMgXII * ne
              + 10**R_MgXIII_to_MgXII_via_e(Tx) * nMgXIII * ne)

dnMgXIII_dt = (
               + 10**R_MgXII_to_MgXIII_via_e(Tx) * nMgXII * ne
               - 10**R_MgXIII_to_MgXII_via_e(Tx) * nMgXIII * ne)

