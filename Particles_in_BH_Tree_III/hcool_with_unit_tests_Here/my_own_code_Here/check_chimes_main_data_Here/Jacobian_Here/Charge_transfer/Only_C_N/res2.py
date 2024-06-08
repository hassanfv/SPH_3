R_HII_to_HI_via_e_caseA = interp1d(Temp, R_HII_to_HI_via_e_caseA_, kind="linear", fill_value="extrapolate") # H CaseA
R_HeII_to_HeI_via_e_caseA = interp1d(Temp, R_HeII_to_HeI_via_e_caseA_, kind="linear", fill_value="extrapolate") # He CaseA
R_HI_to_HII_via_e = interp1d(Temp, R_HI_to_HII_via_e_, kind="linear", fill_value="extrapolate")
R_HI_to_Hm_via_e = interp1d(Temp, R_HI_to_Hm_via_e_, kind="linear", fill_value="extrapolate")
R_Hm_to_HI_via_HI = interp1d(Temp, R_Hm_to_HI_via_HI_, kind="linear", fill_value="extrapolate")
R_Hm_to_HI_via_e = interp1d(Temp, R_Hm_to_HI_via_e_, kind="linear", fill_value="extrapolate")
R_Hm_to_HI_via_HII = interp1d(Temp, R_Hm_to_HI_via_HII_, kind="linear", fill_value="extrapolate")
R_HeI_to_HeII_via_HII = interp1d(Temp, R_HeI_to_HeII_via_HII_, kind="linear", fill_value="extrapolate")
R_HeI_to_HeII_via_e = interp1d(Temp, R_HeI_to_HeII_via_e_, kind="linear", fill_value="extrapolate")
R_HeII_to_HeIII_via_e = interp1d(Temp, R_HeII_to_HeIII_via_e_, kind="linear", fill_value="extrapolate")
R_HeII_to_HeI_via_Hm = interp1d(Temp, R_HeII_to_HeI_via_Hm_, kind="linear", fill_value="extrapolate")
R_HeII_to_HeI_via_HI = interp1d(Temp, R_HeII_to_HeI_via_HI_, kind="linear", fill_value="extrapolate")
R_HeIII_to_HeII_via_HI = interp1d(Temp, R_HeIII_to_HeII_via_HI_, kind="linear", fill_value="extrapolate")
R_HeIII_to_HeII_via_e = interp1d(Temp, R_HeIII_to_HeII_via_e_, kind="linear", fill_value="extrapolate")
