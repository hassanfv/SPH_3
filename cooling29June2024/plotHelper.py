TA = time.time()

print("Running ...")

nH = 1000.0

He_solar = 10**(-1.00)
nHe = He_solar * nH

C_solar = 10**(-3.61)
nC = C_solar * nH

N_solar = 10**(-4.07)
nN = N_solar * nH

O_solar = 10**(-3.31)
nO = O_solar * nH

Ne_solar = 10**(-4.00)
nNe = Ne_solar * nH

Mg_solar = 10**(-4.46)
nMg = Mg_solar * nH

Si_solar = 10**(-4.46)
nSi = Si_solar * nH

S_solar = 10**(-4.73)
nS = S_solar * nH

Ca_solar = 10**(-5.64)
nCa = Ca_solar * nH

Fe_solar = 10**(-4.55)
nFe = Fe_solar * nH

T_i = 10**7.00

nHm_i = 1e-5 * nH
nH0_i = 0.001 * nH
nHp_i = nH - nH0_i

nHe0_i = 0.0001 * nHe
nHep_i = 0.001 * nHe
nHepp_i= nHe - nHe0_i - nHep_i

nC0_i = 1e-6 * nC
nC1_i = 1e-5 * nC
nC2_i = 1e-4 * nC
nC3_i = 1e-3 * nC
nC4_i = 1e-2 * nC
nC5_i = 1e-2 * nC
nCm_i = 1e-6 * nC
nC6_i = nC - (nC0_i + nC1_i + nC2_i + nC3_i + nC4_i + nC5_i + nCm_i)

nN0_i = 1e-7 * nN
nN1_i = 1e-6 * nN
nN2_i = 1e-5 * nN
nN3_i = 1e-4 * nN
nN4_i = 1e-3 * nN
nN5_i = 1e-2 * nN
nN6_i = 1e-2 * nN
nN7_i = nN - (nN0_i + nN1_i + nN2_i + nN3_i + nN4_i + nN5_i + nN6_i)

nO0_i = 1e-8 * nO
nO1_i = 1e-7 * nO
nO2_i = 1e-6 * nO
nO3_i = 1e-5 * nO
nO4_i = 1e-4 * nO
nO5_i = 1e-3 * nO
nO6_i = 1e-2 * nO
nO7_i = 1e-2 * nO
nOm_i = 1e-6 * nO
nO8_i = nO - (nO0_i + nO1_i + nO2_i + nO3_i + nO4_i + nO5_i + nO6_i + nO7_i + nOm_i)

nNe0_i = 1e-10 * nNe
nNe1_i = 1e-9 * nNe
nNe2_i = 1e-8 * nNe
nNe3_i = 1e-7 * nNe
nNe4_i = 1e-6 * nNe
nNe5_i = 1e-5 * nNe
nNe6_i = 1e-4 * nNe
nNe7_i = 1e-3 * nNe
nNe8_i = 1e-2 * nNe
nNe9_i = 1e-2 * nNe
nNe10_i = nNe - (nNe0_i + nNe1_i + nNe2_i + nNe3_i + nNe4_i + nNe5_i + nNe6_i + nNe7_i + nNe8_i + nNe9_i)

nMg0_i = 1e-12 * nMg
nMg1_i = 1e-11 * nMg
nMg2_i = 1e-10 * nMg
nMg3_i = 1e-9 * nMg
nMg4_i = 1e-8 * nMg
nMg5_i = 1e-7 * nMg
nMg6_i = 1e-6 * nMg
nMg7_i = 1e-5 * nMg
nMg8_i = 1e-4 * nMg
nMg9_i = 1e-3 * nMg
nMg10_i = 1e-2 * nMg
nMg11_i = 1e-2 * nMg
nMg12_i = nMg - (nMg0_i + nMg1_i + nMg2_i + nMg3_i + nMg4_i + nMg5_i + nMg6_i + nMg7_i + nMg8_i + nMg9_i + nMg10_i + nMg11_i)

nSi0_i = 1e-14 * nSi
nSi1_i = 1e-13 * nSi
nSi2_i = 1e-12 * nSi
nSi3_i = 1e-11 * nSi
nSi4_i = 1e-10 * nSi
nSi5_i = 1e-9 * nSi
nSi6_i = 1e-8 * nSi
nSi7_i = 1e-7 * nSi
nSi8_i = 1e-6 * nSi
nSi9_i = 1e-5 * nSi
nSi10_i = 1e-4 * nSi
nSi11_i = 1e-3 * nSi
nSi12_i = 1e-2 * nSi
nSi13_i = 1e-2 * nSi
nSi14_i = nSi - (nSi0_i + nSi1_i + nSi2_i + nSi3_i + nSi4_i + nSi5_i + nSi6_i + nSi7_i + nSi8_i + nSi9_i + nSi10_i + nSi11_i + nSi12_i + nSi13_i)

nS0_i = 1e-16 * nS
nS1_i = 1e-15 * nS
nS2_i = 1e-14 * nS
nS3_i = 1e-13 * nS
nS4_i = 1e-12 * nS
nS5_i = 1e-11 * nS
nS6_i = 1e-10 * nS
nS7_i = 1e-9 * nS
nS8_i = 1e-8 * nS
nS9_i = 1e-7 * nS
nS10_i = 1e-6 * nS
nS11_i = 1e-5 * nS
nS12_i = 1e-4 * nS
nS13_i = 1e-3 * nS
nS14_i = 1e-2 * nS
nS15_i = 1e-2 * nS
nS16_i = nS - (nS0_i + nS1_i + nS2_i + nS3_i + nS4_i + nS5_i + nS6_i + nS7_i + nS8_i + nS9_i + nS10_i + nS11_i + nS12_i + nS13_i + nS14_i + nS15_i)

nCa0_i = 1e-20 * nCa
nCa1_i = 1e-19 * nCa
nCa2_i = 1e-18 * nCa
nCa3_i = 1e-17 * nCa
nCa4_i = 1e-16 * nCa
nCa5_i = 1e-15 * nCa
nCa6_i = 1e-14 * nCa
nCa7_i = 1e-13 * nCa
nCa8_i = 1e-12 * nCa
nCa9_i = 1e-11 * nCa
nCa10_i = 1e-10 * nCa
nCa11_i = 1e-9 * nCa
nCa12_i = 1e-8 * nCa
nCa13_i = 1e-7 * nCa
nCa14_i = 1e-6 * nCa
nCa15_i = 1e-5 * nCa
nCa16_i = 1e-4 * nCa
nCa17_i = 1e-3 * nCa
nCa18_i = 1e-2 * nCa
nCa19_i = 1e-2 * nCa
nCa20_i = nCa - (nCa0_i + nCa1_i + nCa2_i + nCa3_i + nCa4_i + nCa5_i + nCa6_i + nCa7_i + nCa8_i + nCa9_i + nCa10_i + nCa11_i + nCa12_i + nCa13_i + nCa14_i + nCa15_i + nCa16_i + nCa17_i + nCa18_i + nCa19_i)

nFe0_i = 1e-26 * nFe
nFe1_i = 1e-25 * nFe
nFe2_i = 1e-24 * nFe
nFe3_i = 1e-23 * nFe
nFe4_i = 1e-22 * nFe
nFe5_i = 1e-21 * nFe
nFe6_i = 1e-20 * nFe
nFe7_i = 1e-19 * nFe
nFe8_i = 1e-18 * nFe
nFe9_i = 1e-17 * nFe
nFe10_i = 1e-16 * nFe
nFe11_i = 1e-15 * nFe
nFe12_i = 1e-14 * nFe
nFe13_i = 1e-13 * nFe
nFe14_i = 1e-12 * nFe
nFe15_i = 1e-11 * nFe
nFe16_i = 1e-10 * nFe
nFe17_i = 1e-9 * nFe
nFe18_i = 1e-8 * nFe
nFe19_i = 1e-7 * nFe
nFe20_i = 1e-6 * nFe
nFe21_i = 1e-5 * nFe
nFe22_i = 1e-4 * nFe
nFe23_i = 1e-3 * nFe
nFe24_i = 1e-2 * nFe
nFe25_i = 1e-2 * nFe
nFe26_i = nFe - (nFe0_i + nFe1_i + nFe2_i + nFe3_i + nFe4_i + nFe5_i + nFe6_i + nFe7_i + nFe8_i + nFe9_i + nFe10_i + nFe11_i + nFe12_i + nFe13_i + nFe14_i + nFe15_i + nFe16_i + nFe17_i + nFe18_i + nFe19_i + nFe20_i + nFe21_i + nFe22_i + nFe23_i + nFe24_i + nFe25_i)

y0 = [
      nH0_i, nHp_i, nHm_i, nHe0_i, nHep_i, nHepp_i, 
      nC0_i, nC1_i, nC2_i, nC3_i, nC4_i, nC5_i, 
      nC6_i, nCm_i, nN0_i, nN1_i, nN2_i, nN3_i, 
      nN4_i, nN5_i, nN6_i, nN7_i, nO0_i, nO1_i, 
      nO2_i, nO3_i, nO4_i, nO5_i, nO6_i, nO7_i, 
      nO8_i, nOm_i, nNe0_i, nNe1_i, nNe2_i, nNe3_i, 
      nNe4_i, nNe5_i, nNe6_i, nNe7_i, nNe8_i, nNe9_i, 
      nNe10_i, nMg0_i, nMg1_i, nMg2_i, nMg3_i, nMg4_i, 
      nMg5_i, nMg6_i, nMg7_i, nMg8_i, nMg9_i, nMg10_i, 
      nMg11_i, nMg12_i, nSi0_i, nSi1_i, nSi2_i, nSi3_i, 
      nSi4_i, nSi5_i, nSi6_i, nSi7_i, nSi8_i, nSi9_i, 
      nSi10_i, nSi11_i, nSi12_i, nSi13_i, nSi14_i, nS0_i, 
      nS1_i, nS2_i, nS3_i, nS4_i, nS5_i, nS6_i, 
      nS7_i, nS8_i, nS9_i, nS10_i, nS11_i, nS12_i, 
      nS13_i, nS14_i, nS15_i, nS16_i, nCa0_i, nCa1_i, 
      nCa2_i, nCa3_i, nCa4_i, nCa5_i, nCa6_i, nCa7_i, 
      nCa8_i, nCa9_i, nCa10_i, nCa11_i, nCa12_i, nCa13_i, 
      nCa14_i, nCa15_i, nCa16_i, nCa17_i, nCa18_i, nCa19_i, 
      nCa20_i, nFe0_i, nFe1_i, nFe2_i, nFe3_i, nFe4_i, 
      nFe5_i, nFe6_i, nFe7_i, nFe8_i, nFe9_i, nFe10_i, 
      nFe11_i, nFe12_i, nFe13_i, nFe14_i, nFe15_i, nFe16_i, 
      nFe17_i, nFe18_i, nFe19_i, nFe20_i, nFe21_i, nFe22_i, 
      nFe23_i, nFe24_i, nFe25_i, nFe26_i, 
      T_i
      ]

A_v = 1.0
G0 = 0.01
dust_ratio = 0.01

t_span = (1*3.16e7, 10000*3.16e7)

solution = solve_ivp(func, t_span, y0, method="LSODA", dense_output=True)

t = np.linspace(t_span[0], t_span[1], 50000) # This 10000 is not years, it is the number of points in linspace !!!!
y = solution.sol(t)

t_yrs = t / 3.16e7

nH0  = y[0, :]
nHp  = y[1, :]
nHm  = y[2, :]
nHe0 = y[3, :]
nHep = y[4, :]
nHepp= y[5, :]

nC0 = y[6, :]
nC1 = y[7, :]
nC2 = y[8, :]
nC3 = y[9, :]
nC4 = y[10, :]
nC5 = y[11, :]
nC6 = y[12, :]
nCm = y[13, :]


nN0 = y[14, :]
nN1 = y[15, :]
nN2 = y[16, :]
nN3 = y[17, :]
nN4 = y[18, :]
nN5 = y[19, :]
nN6 = y[20, :]
nN7 = y[21, :]


nO0 = y[22, :]
nO1 = y[23, :]
nO2 = y[24, :]
nO3 = y[25, :]
nO4 = y[26, :]
nO5 = y[27, :]
nO6 = y[28, :]
nO7 = y[29, :]
nO8 = y[30, :]
nOm = y[31, :]


nNe0 = y[32, :]
nNe1 = y[33, :]
nNe2 = y[34, :]
nNe3 = y[35, :]
nNe4 = y[36, :]
nNe5 = y[37, :]
nNe6 = y[38, :]
nNe7 = y[39, :]
nNe8 = y[40, :]
nNe9 = y[41, :]
nNe10 = y[42, :]


nMg0 = y[43, :]
nMg1 = y[44, :]
nMg2 = y[45, :]
nMg3 = y[46, :]
nMg4 = y[47, :]
nMg5 = y[48, :]
nMg6 = y[49, :]
nMg7 = y[50, :]
nMg8 = y[51, :]
nMg9 = y[52, :]
nMg10 = y[53, :]
nMg11 = y[54, :]
nMg12 = y[55, :]


nSi0 = y[56, :]
nSi1 = y[57, :]
nSi2 = y[58, :]
nSi3 = y[59, :]
nSi4 = y[60, :]
nSi5 = y[61, :]
nSi6 = y[62, :]
nSi7 = y[63, :]
nSi8 = y[64, :]
nSi9 = y[65, :]
nSi10 = y[66, :]
nSi11 = y[67, :]
nSi12 = y[68, :]
nSi13 = y[69, :]
nSi14 = y[70, :]


nS0 = y[71, :]
nS1 = y[72, :]
nS2 = y[73, :]
nS3 = y[74, :]
nS4 = y[75, :]
nS5 = y[76, :]
nS6 = y[77, :]
nS7 = y[78, :]
nS8 = y[79, :]
nS9 = y[80, :]
nS10 = y[81, :]
nS11 = y[82, :]
nS12 = y[83, :]
nS13 = y[84, :]
nS14 = y[85, :]
nS15 = y[86, :]
nS16 = y[87, :]


nCa0 = y[88, :]
nCa1 = y[89, :]
nCa2 = y[90, :]
nCa3 = y[91, :]
nCa4 = y[92, :]
nCa5 = y[93, :]
nCa6 = y[94, :]
nCa7 = y[95, :]
nCa8 = y[96, :]
nCa9 = y[97, :]
nCa10 = y[98, :]
nCa11 = y[99, :]
nCa12 = y[100, :]
nCa13 = y[101, :]
nCa14 = y[102, :]
nCa15 = y[103, :]
nCa16 = y[104, :]
nCa17 = y[105, :]
nCa18 = y[106, :]
nCa19 = y[107, :]
nCa20 = y[108, :]


nFe0 = y[109, :]
nFe1 = y[110, :]
nFe2 = y[111, :]
nFe3 = y[112, :]
nFe4 = y[113, :]
nFe5 = y[114, :]
nFe6 = y[115, :]
nFe7 = y[116, :]
nFe8 = y[117, :]
nFe9 = y[118, :]
nFe10 = y[119, :]
nFe11 = y[120, :]
nFe12 = y[121, :]
nFe13 = y[122, :]
nFe14 = y[123, :]
nFe15 = y[124, :]
nFe16 = y[125, :]
nFe17 = y[126, :]
nFe18 = y[127, :]
nFe19 = y[128, :]
nFe20 = y[129, :]
nFe21 = y[130, :]
nFe22 = y[131, :]
nFe23 = y[132, :]
nFe24 = y[133, :]
nFe25 = y[134, :]
nFe26 = y[135, :]


T = y[136, :]

with open ("chimesData.pkl", "rb") as f:
  data = pickle.load(f)

TEvolx = data["TempEvol"]
AbundEvol = data["chimesAbundEvol"]
t_Arr_in_yrsx = data["t_Arr_in_yrs"]

nH0x   = AbundEvol[1, :]
nHpx   = AbundEvol[2, :]
nHe0x  = AbundEvol[4, :]
nHepx  = AbundEvol[5, :]
nHeppx = AbundEvol[6, :]

nC0x = AbundEvol[7, :]
nC1x = AbundEvol[8, :]
nC2x = AbundEvol[9, :]
nC3x = AbundEvol[10, :]
nC4x = AbundEvol[11, :]
nC5x = AbundEvol[12, :]
nC6x = AbundEvol[13, :]
nCmx = AbundEvol[14, :]
nCx = nC0x + nC1x + nC2x + nC3x + nC4x + nC5x + nC6x + nCmx 

nN0x = AbundEvol[15, :]
nN1x = AbundEvol[16, :]
nN2x = AbundEvol[17, :]
nN3x = AbundEvol[18, :]
nN4x = AbundEvol[19, :]
nN5x = AbundEvol[20, :]
nN6x = AbundEvol[21, :]
nN7x = AbundEvol[22, :]
nNx = nN0x + nN1x + nN2x + nN3x + nN4x + nN5x + nN6x + nN7x 

nO0x = AbundEvol[23, :]
nO1x = AbundEvol[24, :]
nO2x = AbundEvol[25, :]
nO3x = AbundEvol[26, :]
nO4x = AbundEvol[27, :]
nO5x = AbundEvol[28, :]
nO6x = AbundEvol[29, :]
nO7x = AbundEvol[30, :]
nO8x = AbundEvol[31, :]
nOmx = AbundEvol[32, :]
nOx = nO0x + nO1x + nO2x + nO3x + nO4x + nO5x + nO6x + nO7x + nO8x + nOmx 

nNe0x = AbundEvol[33, :]
nNe1x = AbundEvol[34, :]
nNe2x = AbundEvol[35, :]
nNe3x = AbundEvol[36, :]
nNe4x = AbundEvol[37, :]
nNe5x = AbundEvol[38, :]
nNe6x = AbundEvol[39, :]
nNe7x = AbundEvol[40, :]
nNe8x = AbundEvol[41, :]
nNe9x = AbundEvol[42, :]
nNe10x = AbundEvol[43, :]
nNex = nNe0x + nNe1x + nNe2x + nNe3x + nNe4x + nNe5x + nNe6x + nNe7x + nNe8x + nNe9x + nNe10x 

nMg0x = AbundEvol[44, :]
nMg1x = AbundEvol[45, :]
nMg2x = AbundEvol[46, :]
nMg3x = AbundEvol[47, :]
nMg4x = AbundEvol[48, :]
nMg5x = AbundEvol[49, :]
nMg6x = AbundEvol[50, :]
nMg7x = AbundEvol[51, :]
nMg8x = AbundEvol[52, :]
nMg9x = AbundEvol[53, :]
nMg10x = AbundEvol[54, :]
nMg11x = AbundEvol[55, :]
nMg12x = AbundEvol[56, :]
nMgx = nMg0x + nMg1x + nMg2x + nMg3x + nMg4x + nMg5x + nMg6x + nMg7x + nMg8x + nMg9x + nMg10x + nMg11x + nMg12x 

nSi0x = AbundEvol[57, :]
nSi1x = AbundEvol[58, :]
nSi2x = AbundEvol[59, :]
nSi3x = AbundEvol[60, :]
nSi4x = AbundEvol[61, :]
nSi5x = AbundEvol[62, :]
nSi6x = AbundEvol[63, :]
nSi7x = AbundEvol[64, :]
nSi8x = AbundEvol[65, :]
nSi9x = AbundEvol[66, :]
nSi10x = AbundEvol[67, :]
nSi11x = AbundEvol[68, :]
nSi12x = AbundEvol[69, :]
nSi13x = AbundEvol[70, :]
nSi14x = AbundEvol[71, :]
nSix = nSi0x + nSi1x + nSi2x + nSi3x + nSi4x + nSi5x + nSi6x + nSi7x + nSi8x + nSi9x + nSi10x + nSi11x + nSi12x + nSi13x + nSi14x 

nS0x = AbundEvol[72, :]
nS1x = AbundEvol[73, :]
nS2x = AbundEvol[74, :]
nS3x = AbundEvol[75, :]
nS4x = AbundEvol[76, :]
nS5x = AbundEvol[77, :]
nS6x = AbundEvol[78, :]
nS7x = AbundEvol[79, :]
nS8x = AbundEvol[80, :]
nS9x = AbundEvol[81, :]
nS10x = AbundEvol[82, :]
nS11x = AbundEvol[83, :]
nS12x = AbundEvol[84, :]
nS13x = AbundEvol[85, :]
nS14x = AbundEvol[86, :]
nS15x = AbundEvol[87, :]
nS16x = AbundEvol[88, :]
nSx = nS0x + nS1x + nS2x + nS3x + nS4x + nS5x + nS6x + nS7x + nS8x + nS9x + nS10x + nS11x + nS12x + nS13x + nS14x + nS15x + nS16x 

nCa0x = AbundEvol[89, :]
nCa1x = AbundEvol[90, :]
nCa2x = AbundEvol[91, :]
nCa3x = AbundEvol[92, :]
nCa4x = AbundEvol[93, :]
nCa5x = AbundEvol[94, :]
nCa6x = AbundEvol[95, :]
nCa7x = AbundEvol[96, :]
nCa8x = AbundEvol[97, :]
nCa9x = AbundEvol[98, :]
nCa10x = AbundEvol[99, :]
nCa11x = AbundEvol[100, :]
nCa12x = AbundEvol[101, :]
nCa13x = AbundEvol[102, :]
nCa14x = AbundEvol[103, :]
nCa15x = AbundEvol[104, :]
nCa16x = AbundEvol[105, :]
nCa17x = AbundEvol[106, :]
nCa18x = AbundEvol[107, :]
nCa19x = AbundEvol[108, :]
nCa20x = AbundEvol[109, :]
nCax = nCa0x + nCa1x + nCa2x + nCa3x + nCa4x + nCa5x + nCa6x + nCa7x + nCa8x + nCa9x + nCa10x + nCa11x + nCa12x + nCa13x + nCa14x + nCa15x + nCa16x + nCa17x + nCa18x + nCa19x + nCa20x 

nFe0x = AbundEvol[111, :]
nFe1x = AbundEvol[112, :]
nFe2x = AbundEvol[113, :]
nFe3x = AbundEvol[114, :]
nFe4x = AbundEvol[115, :]
nFe5x = AbundEvol[116, :]
nFe6x = AbundEvol[117, :]
nFe7x = AbundEvol[118, :]
nFe8x = AbundEvol[119, :]
nFe9x = AbundEvol[120, :]
nFe10x = AbundEvol[121, :]
nFe11x = AbundEvol[122, :]
nFe12x = AbundEvol[123, :]
nFe13x = AbundEvol[124, :]
nFe14x = AbundEvol[125, :]
nFe15x = AbundEvol[126, :]
nFe16x = AbundEvol[127, :]
nFe17x = AbundEvol[128, :]
nFe18x = AbundEvol[129, :]
nFe19x = AbundEvol[130, :]
nFe20x = AbundEvol[131, :]
nFe21x = AbundEvol[132, :]
nFe22x = AbundEvol[133, :]
nFe23x = AbundEvol[134, :]
nFe24x = AbundEvol[135, :]
nFe25x = AbundEvol[136, :]
nFe26x = AbundEvol[137, :]
nFex = nFe0x + nFe1x + nFe2x + nFe3x + nFe4x + nFe5x + nFe6x + nFe7x + nFe8x + nFe9x + nFe10x + nFe11x + nFe12x + nFe13x + nFe14x + nFe15x + nFe16x + nFe17x + nFe18x + nFe19x + nFe20x + nFe21x + nFe22x + nFe23x + nFe24x + nFe25x + nFe26x 

plt.figure(figsize=(10, 5))
plt.scatter(t_yrs, np.log10(T), s=2, color="k", label="my code")
plt.scatter(t_Arr_in_yrsx, np.log10(TEvolx), s=2, color="orange", label="chimes result", linestyle="--")
plt.xlim(0, 10000)
plt.ylim(1, 8)
plt.legend()
plt.savefig("Temp_vs_time.png", dpi = 300)
plt.close()

plt.plot(t_yrs, nHe0, color = 'r', label = 'nHe0')
plt.plot(t_yrs, nHep, color = 'g', label = 'nHep')
plt.plot(t_yrs, nHepp, color = 'b', label = 'nHepp')
plt.plot(t_Arr_in_yrsx, nHe0x, color = 'r', label = 'nHe0 - chimes', linestyle = ':')
plt.plot(t_Arr_in_yrsx, nHepx, color = 'g', label = 'nHep - chimes', linestyle = ':')
plt.plot(t_Arr_in_yrsx, nHeppx, color = 'b', label = 'nHepp - chimes', linestyle = ':')
plt.xlim(0, 10000)
plt.ylim(1e-8, 300)
plt.yscale('log')
plt.title('solve_ivp')
plt.legend()
plt.savefig("nH_He_vs_time.png", dpi = 300)
plt.close()

plt.plot(T, nC0/nC, label = 'nC0', color = '#1f77b4')
plt.plot(TEvolx, nC0x/nCx, label = 'nC0x', color = '#1f77b4', linestyle = ':')
plt.plot(T, nC1/nC, label = 'nC1', color = '#ff7f0e')
plt.plot(TEvolx, nC1x/nCx, label = 'nC1x', color = '#ff7f0e', linestyle = ':')
plt.plot(T, nC2/nC, label = 'nC2', color = '#2ca02c')
plt.plot(TEvolx, nC2x/nCx, label = 'nC2x', color = '#2ca02c', linestyle = ':')
plt.plot(T, nC3/nC, label = 'nC3', color = '#d62728')
plt.plot(TEvolx, nC3x/nCx, label = 'nC3x', color = '#d62728', linestyle = ':')
plt.plot(T, nC4/nC, label = 'nC4', color = '#9467bd')
plt.plot(TEvolx, nC4x/nCx, label = 'nC4x', color = '#9467bd', linestyle = ':')
plt.plot(T, nC5/nC, label = 'nC5', color = '#8c564b')
plt.plot(TEvolx, nC5x/nCx, label = 'nC5x', color = '#8c564b', linestyle = ':')
plt.plot(T, nC6/nC, label = 'nC6', color = '#e377c2')
plt.plot(TEvolx, nC6x/nCx, label = 'nC6x', color = '#e377c2', linestyle = ':')
plt.plot(T, nCm/nC, label = 'nCm', color = '#7f7f7f')
plt.plot(TEvolx, nCmx/nCx, label = 'nCmx', color = '#7f7f7f', linestyle = ':')
plt.yscale('log')
plt.xscale('log')
plt.ylim(2e-3, 1.2)
plt.xlim(1e4, 1e7)
plt.savefig("nC_vs_time.png", dpi = 300)
plt.close()

plt.plot(T, nN0/nN, label = 'nN0', color = '#1f77b4')
plt.plot(TEvolx, nN0x/nNx, label = 'nN0x', color = '#1f77b4', linestyle = ':')
plt.plot(T, nN1/nN, label = 'nN1', color = '#ff7f0e')
plt.plot(TEvolx, nN1x/nNx, label = 'nN1x', color = '#ff7f0e', linestyle = ':')
plt.plot(T, nN2/nN, label = 'nN2', color = '#2ca02c')
plt.plot(TEvolx, nN2x/nNx, label = 'nN2x', color = '#2ca02c', linestyle = ':')
plt.plot(T, nN3/nN, label = 'nN3', color = '#d62728')
plt.plot(TEvolx, nN3x/nNx, label = 'nN3x', color = '#d62728', linestyle = ':')
plt.plot(T, nN4/nN, label = 'nN4', color = '#9467bd')
plt.plot(TEvolx, nN4x/nNx, label = 'nN4x', color = '#9467bd', linestyle = ':')
plt.plot(T, nN5/nN, label = 'nN5', color = '#8c564b')
plt.plot(TEvolx, nN5x/nNx, label = 'nN5x', color = '#8c564b', linestyle = ':')
plt.plot(T, nN6/nN, label = 'nN6', color = '#e377c2')
plt.plot(TEvolx, nN6x/nNx, label = 'nN6x', color = '#e377c2', linestyle = ':')
plt.plot(T, nN7/nN, label = 'nN7', color = '#7f7f7f')
plt.plot(TEvolx, nN7x/nNx, label = 'nN7x', color = '#7f7f7f', linestyle = ':')
plt.yscale('log')
plt.xscale('log')
plt.ylim(2e-3, 1.2)
plt.xlim(1e4, 1e7)
plt.savefig("nN_vs_time.png", dpi = 300)
plt.close()

plt.plot(T, nO0/nO, label = 'nO0', color = '#1f77b4')
plt.plot(TEvolx, nO0x/nOx, label = 'nO0x', color = '#1f77b4', linestyle = ':')
plt.plot(T, nO1/nO, label = 'nO1', color = '#ff7f0e')
plt.plot(TEvolx, nO1x/nOx, label = 'nO1x', color = '#ff7f0e', linestyle = ':')
plt.plot(T, nO2/nO, label = 'nO2', color = '#2ca02c')
plt.plot(TEvolx, nO2x/nOx, label = 'nO2x', color = '#2ca02c', linestyle = ':')
plt.plot(T, nO3/nO, label = 'nO3', color = '#d62728')
plt.plot(TEvolx, nO3x/nOx, label = 'nO3x', color = '#d62728', linestyle = ':')
plt.plot(T, nO4/nO, label = 'nO4', color = '#9467bd')
plt.plot(TEvolx, nO4x/nOx, label = 'nO4x', color = '#9467bd', linestyle = ':')
plt.plot(T, nO5/nO, label = 'nO5', color = '#8c564b')
plt.plot(TEvolx, nO5x/nOx, label = 'nO5x', color = '#8c564b', linestyle = ':')
plt.plot(T, nO6/nO, label = 'nO6', color = '#e377c2')
plt.plot(TEvolx, nO6x/nOx, label = 'nO6x', color = '#e377c2', linestyle = ':')
plt.plot(T, nO7/nO, label = 'nO7', color = '#7f7f7f')
plt.plot(TEvolx, nO7x/nOx, label = 'nO7x', color = '#7f7f7f', linestyle = ':')
plt.plot(T, nO8/nO, label = 'nO8', color = '#bcbd22')
plt.plot(TEvolx, nO8x/nOx, label = 'nO8x', color = '#bcbd22', linestyle = ':')
plt.plot(T, nOm/nO, label = 'nOm', color = '#17becf')
plt.plot(TEvolx, nOmx/nOx, label = 'nOmx', color = '#17becf', linestyle = ':')
plt.yscale('log')
plt.xscale('log')
plt.ylim(2e-3, 1.2)
plt.xlim(1e4, 1e7)
plt.savefig("nO_vs_time.png", dpi = 300)
plt.close()

plt.plot(T, nNe0/nNe, label = 'nNe0', color = '#1f77b4')
plt.plot(TEvolx, nNe0x/nNex, label = 'nNe0x', color = '#1f77b4', linestyle = ':')
plt.plot(T, nNe1/nNe, label = 'nNe1', color = '#ff7f0e')
plt.plot(TEvolx, nNe1x/nNex, label = 'nNe1x', color = '#ff7f0e', linestyle = ':')
plt.plot(T, nNe2/nNe, label = 'nNe2', color = '#2ca02c')
plt.plot(TEvolx, nNe2x/nNex, label = 'nNe2x', color = '#2ca02c', linestyle = ':')
plt.plot(T, nNe3/nNe, label = 'nNe3', color = '#d62728')
plt.plot(TEvolx, nNe3x/nNex, label = 'nNe3x', color = '#d62728', linestyle = ':')
plt.plot(T, nNe4/nNe, label = 'nNe4', color = '#9467bd')
plt.plot(TEvolx, nNe4x/nNex, label = 'nNe4x', color = '#9467bd', linestyle = ':')
plt.plot(T, nNe5/nNe, label = 'nNe5', color = '#8c564b')
plt.plot(TEvolx, nNe5x/nNex, label = 'nNe5x', color = '#8c564b', linestyle = ':')
plt.plot(T, nNe6/nNe, label = 'nNe6', color = '#e377c2')
plt.plot(TEvolx, nNe6x/nNex, label = 'nNe6x', color = '#e377c2', linestyle = ':')
plt.plot(T, nNe7/nNe, label = 'nNe7', color = '#7f7f7f')
plt.plot(TEvolx, nNe7x/nNex, label = 'nNe7x', color = '#7f7f7f', linestyle = ':')
plt.plot(T, nNe8/nNe, label = 'nNe8', color = '#bcbd22')
plt.plot(TEvolx, nNe8x/nNex, label = 'nNe8x', color = '#bcbd22', linestyle = ':')
plt.plot(T, nNe9/nNe, label = 'nNe9', color = '#17becf')
plt.plot(TEvolx, nNe9x/nNex, label = 'nNe9x', color = '#17becf', linestyle = ':')
plt.plot(T, nNe10/nNe, label = 'nNe10', color = '#1f77b4')
plt.plot(TEvolx, nNe10x/nNex, label = 'nNe10x', color = '#1f77b4', linestyle = ':')
plt.yscale('log')
plt.xscale('log')
plt.ylim(2e-3, 1.2)
plt.xlim(1e4, 1e7)
plt.savefig("nNe_vs_time.png", dpi = 300)
plt.close()

plt.plot(T, nMg0/nMg, label = 'nMg0', color = '#1f77b4')
plt.plot(TEvolx, nMg0x/nMgx, label = 'nMg0x', color = '#1f77b4', linestyle = ':')
plt.plot(T, nMg1/nMg, label = 'nMg1', color = '#ff7f0e')
plt.plot(TEvolx, nMg1x/nMgx, label = 'nMg1x', color = '#ff7f0e', linestyle = ':')
plt.plot(T, nMg2/nMg, label = 'nMg2', color = '#2ca02c')
plt.plot(TEvolx, nMg2x/nMgx, label = 'nMg2x', color = '#2ca02c', linestyle = ':')
plt.plot(T, nMg3/nMg, label = 'nMg3', color = '#d62728')
plt.plot(TEvolx, nMg3x/nMgx, label = 'nMg3x', color = '#d62728', linestyle = ':')
plt.plot(T, nMg4/nMg, label = 'nMg4', color = '#9467bd')
plt.plot(TEvolx, nMg4x/nMgx, label = 'nMg4x', color = '#9467bd', linestyle = ':')
plt.plot(T, nMg5/nMg, label = 'nMg5', color = '#8c564b')
plt.plot(TEvolx, nMg5x/nMgx, label = 'nMg5x', color = '#8c564b', linestyle = ':')
plt.plot(T, nMg6/nMg, label = 'nMg6', color = '#e377c2')
plt.plot(TEvolx, nMg6x/nMgx, label = 'nMg6x', color = '#e377c2', linestyle = ':')
plt.plot(T, nMg7/nMg, label = 'nMg7', color = '#7f7f7f')
plt.plot(TEvolx, nMg7x/nMgx, label = 'nMg7x', color = '#7f7f7f', linestyle = ':')
plt.plot(T, nMg8/nMg, label = 'nMg8', color = '#bcbd22')
plt.plot(TEvolx, nMg8x/nMgx, label = 'nMg8x', color = '#bcbd22', linestyle = ':')
plt.plot(T, nMg9/nMg, label = 'nMg9', color = '#17becf')
plt.plot(TEvolx, nMg9x/nMgx, label = 'nMg9x', color = '#17becf', linestyle = ':')
plt.plot(T, nMg10/nMg, label = 'nMg10', color = '#1f77b4')
plt.plot(TEvolx, nMg10x/nMgx, label = 'nMg10x', color = '#1f77b4', linestyle = ':')
plt.plot(T, nMg11/nMg, label = 'nMg11', color = '#aec7e8')
plt.plot(TEvolx, nMg11x/nMgx, label = 'nMg11x', color = '#aec7e8', linestyle = ':')
plt.plot(T, nMg12/nMg, label = 'nMg12', color = '#ffbb78')
plt.plot(TEvolx, nMg12x/nMgx, label = 'nMg12x', color = '#ffbb78', linestyle = ':')
plt.yscale('log')
plt.xscale('log')
plt.ylim(2e-3, 1.2)
plt.xlim(1e4, 1e7)
plt.savefig("nMg_vs_time.png", dpi = 300)
plt.close()

plt.plot(T, nSi0/nSi, label = 'nSi0', color = '#1f77b4')
plt.plot(TEvolx, nSi0x/nSix, label = 'nSi0x', color = '#1f77b4', linestyle = ':')
plt.plot(T, nSi1/nSi, label = 'nSi1', color = '#ff7f0e')
plt.plot(TEvolx, nSi1x/nSix, label = 'nSi1x', color = '#ff7f0e', linestyle = ':')
plt.plot(T, nSi2/nSi, label = 'nSi2', color = '#2ca02c')
plt.plot(TEvolx, nSi2x/nSix, label = 'nSi2x', color = '#2ca02c', linestyle = ':')
plt.plot(T, nSi3/nSi, label = 'nSi3', color = '#d62728')
plt.plot(TEvolx, nSi3x/nSix, label = 'nSi3x', color = '#d62728', linestyle = ':')
plt.plot(T, nSi4/nSi, label = 'nSi4', color = '#9467bd')
plt.plot(TEvolx, nSi4x/nSix, label = 'nSi4x', color = '#9467bd', linestyle = ':')
plt.plot(T, nSi5/nSi, label = 'nSi5', color = '#8c564b')
plt.plot(TEvolx, nSi5x/nSix, label = 'nSi5x', color = '#8c564b', linestyle = ':')
plt.plot(T, nSi6/nSi, label = 'nSi6', color = '#e377c2')
plt.plot(TEvolx, nSi6x/nSix, label = 'nSi6x', color = '#e377c2', linestyle = ':')
plt.plot(T, nSi7/nSi, label = 'nSi7', color = '#7f7f7f')
plt.plot(TEvolx, nSi7x/nSix, label = 'nSi7x', color = '#7f7f7f', linestyle = ':')
plt.plot(T, nSi8/nSi, label = 'nSi8', color = '#bcbd22')
plt.plot(TEvolx, nSi8x/nSix, label = 'nSi8x', color = '#bcbd22', linestyle = ':')
plt.plot(T, nSi9/nSi, label = 'nSi9', color = '#17becf')
plt.plot(TEvolx, nSi9x/nSix, label = 'nSi9x', color = '#17becf', linestyle = ':')
plt.plot(T, nSi10/nSi, label = 'nSi10', color = '#1f77b4')
plt.plot(TEvolx, nSi10x/nSix, label = 'nSi10x', color = '#1f77b4', linestyle = ':')
plt.plot(T, nSi11/nSi, label = 'nSi11', color = '#aec7e8')
plt.plot(TEvolx, nSi11x/nSix, label = 'nSi11x', color = '#aec7e8', linestyle = ':')
plt.plot(T, nSi12/nSi, label = 'nSi12', color = '#ffbb78')
plt.plot(TEvolx, nSi12x/nSix, label = 'nSi12x', color = '#ffbb78', linestyle = ':')
plt.plot(T, nSi13/nSi, label = 'nSi13', color = '#98df8a')
plt.plot(TEvolx, nSi13x/nSix, label = 'nSi13x', color = '#98df8a', linestyle = ':')
plt.plot(T, nSi14/nSi, label = 'nSi14', color = '#ff9896')
plt.plot(TEvolx, nSi14x/nSix, label = 'nSi14x', color = '#ff9896', linestyle = ':')
plt.yscale('log')
plt.xscale('log')
plt.ylim(2e-3, 1.2)
plt.xlim(1e4, 1e7)
plt.savefig("nSi_vs_time.png", dpi = 300)
plt.close()

plt.plot(T, nS0/nS, label = 'nS0', color = '#1f77b4')
plt.plot(TEvolx, nS0x/nSx, label = 'nS0x', color = '#1f77b4', linestyle = ':')
plt.plot(T, nS1/nS, label = 'nS1', color = '#ff7f0e')
plt.plot(TEvolx, nS1x/nSx, label = 'nS1x', color = '#ff7f0e', linestyle = ':')
plt.plot(T, nS2/nS, label = 'nS2', color = '#2ca02c')
plt.plot(TEvolx, nS2x/nSx, label = 'nS2x', color = '#2ca02c', linestyle = ':')
plt.plot(T, nS3/nS, label = 'nS3', color = '#d62728')
plt.plot(TEvolx, nS3x/nSx, label = 'nS3x', color = '#d62728', linestyle = ':')
plt.plot(T, nS4/nS, label = 'nS4', color = '#9467bd')
plt.plot(TEvolx, nS4x/nSx, label = 'nS4x', color = '#9467bd', linestyle = ':')
plt.plot(T, nS5/nS, label = 'nS5', color = '#8c564b')
plt.plot(TEvolx, nS5x/nSx, label = 'nS5x', color = '#8c564b', linestyle = ':')
plt.plot(T, nS6/nS, label = 'nS6', color = '#e377c2')
plt.plot(TEvolx, nS6x/nSx, label = 'nS6x', color = '#e377c2', linestyle = ':')
plt.plot(T, nS7/nS, label = 'nS7', color = '#7f7f7f')
plt.plot(TEvolx, nS7x/nSx, label = 'nS7x', color = '#7f7f7f', linestyle = ':')
plt.plot(T, nS8/nS, label = 'nS8', color = '#bcbd22')
plt.plot(TEvolx, nS8x/nSx, label = 'nS8x', color = '#bcbd22', linestyle = ':')
plt.plot(T, nS9/nS, label = 'nS9', color = '#17becf')
plt.plot(TEvolx, nS9x/nSx, label = 'nS9x', color = '#17becf', linestyle = ':')
plt.plot(T, nS10/nS, label = 'nS10', color = '#1f77b4')
plt.plot(TEvolx, nS10x/nSx, label = 'nS10x', color = '#1f77b4', linestyle = ':')
plt.plot(T, nS11/nS, label = 'nS11', color = '#aec7e8')
plt.plot(TEvolx, nS11x/nSx, label = 'nS11x', color = '#aec7e8', linestyle = ':')
plt.plot(T, nS12/nS, label = 'nS12', color = '#ffbb78')
plt.plot(TEvolx, nS12x/nSx, label = 'nS12x', color = '#ffbb78', linestyle = ':')
plt.plot(T, nS13/nS, label = 'nS13', color = '#98df8a')
plt.plot(TEvolx, nS13x/nSx, label = 'nS13x', color = '#98df8a', linestyle = ':')
plt.plot(T, nS14/nS, label = 'nS14', color = '#ff9896')
plt.plot(TEvolx, nS14x/nSx, label = 'nS14x', color = '#ff9896', linestyle = ':')
plt.plot(T, nS15/nS, label = 'nS15', color = '#c5b0d5')
plt.plot(TEvolx, nS15x/nSx, label = 'nS15x', color = '#c5b0d5', linestyle = ':')
plt.plot(T, nS16/nS, label = 'nS16', color = '#c49c94')
plt.plot(TEvolx, nS16x/nSx, label = 'nS16x', color = '#c49c94', linestyle = ':')
plt.yscale('log')
plt.xscale('log')
plt.ylim(2e-3, 1.2)
plt.xlim(1e4, 1e7)
plt.savefig("nS_vs_time.png", dpi = 300)
plt.close()

plt.plot(T, nCa0/nCa, label = 'nCa0', color = '#1f77b4')
plt.plot(TEvolx, nCa0x/nCax, label = 'nCa0x', color = '#1f77b4', linestyle = ':')
plt.plot(T, nCa1/nCa, label = 'nCa1', color = '#ff7f0e')
plt.plot(TEvolx, nCa1x/nCax, label = 'nCa1x', color = '#ff7f0e', linestyle = ':')
plt.plot(T, nCa2/nCa, label = 'nCa2', color = '#2ca02c')
plt.plot(TEvolx, nCa2x/nCax, label = 'nCa2x', color = '#2ca02c', linestyle = ':')
plt.plot(T, nCa3/nCa, label = 'nCa3', color = '#d62728')
plt.plot(TEvolx, nCa3x/nCax, label = 'nCa3x', color = '#d62728', linestyle = ':')
plt.plot(T, nCa4/nCa, label = 'nCa4', color = '#9467bd')
plt.plot(TEvolx, nCa4x/nCax, label = 'nCa4x', color = '#9467bd', linestyle = ':')
plt.plot(T, nCa5/nCa, label = 'nCa5', color = '#8c564b')
plt.plot(TEvolx, nCa5x/nCax, label = 'nCa5x', color = '#8c564b', linestyle = ':')
plt.plot(T, nCa6/nCa, label = 'nCa6', color = '#e377c2')
plt.plot(TEvolx, nCa6x/nCax, label = 'nCa6x', color = '#e377c2', linestyle = ':')
plt.plot(T, nCa7/nCa, label = 'nCa7', color = '#7f7f7f')
plt.plot(TEvolx, nCa7x/nCax, label = 'nCa7x', color = '#7f7f7f', linestyle = ':')
plt.plot(T, nCa8/nCa, label = 'nCa8', color = '#bcbd22')
plt.plot(TEvolx, nCa8x/nCax, label = 'nCa8x', color = '#bcbd22', linestyle = ':')
plt.plot(T, nCa9/nCa, label = 'nCa9', color = '#17becf')
plt.plot(TEvolx, nCa9x/nCax, label = 'nCa9x', color = '#17becf', linestyle = ':')
plt.plot(T, nCa10/nCa, label = 'nCa10', color = '#1f77b4')
plt.plot(TEvolx, nCa10x/nCax, label = 'nCa10x', color = '#1f77b4', linestyle = ':')
plt.plot(T, nCa11/nCa, label = 'nCa11', color = '#aec7e8')
plt.plot(TEvolx, nCa11x/nCax, label = 'nCa11x', color = '#aec7e8', linestyle = ':')
plt.plot(T, nCa12/nCa, label = 'nCa12', color = '#ffbb78')
plt.plot(TEvolx, nCa12x/nCax, label = 'nCa12x', color = '#ffbb78', linestyle = ':')
plt.plot(T, nCa13/nCa, label = 'nCa13', color = '#98df8a')
plt.plot(TEvolx, nCa13x/nCax, label = 'nCa13x', color = '#98df8a', linestyle = ':')
plt.plot(T, nCa14/nCa, label = 'nCa14', color = '#ff9896')
plt.plot(TEvolx, nCa14x/nCax, label = 'nCa14x', color = '#ff9896', linestyle = ':')
plt.plot(T, nCa15/nCa, label = 'nCa15', color = '#c5b0d5')
plt.plot(TEvolx, nCa15x/nCax, label = 'nCa15x', color = '#c5b0d5', linestyle = ':')
plt.plot(T, nCa16/nCa, label = 'nCa16', color = '#c49c94')
plt.plot(TEvolx, nCa16x/nCax, label = 'nCa16x', color = '#c49c94', linestyle = ':')
plt.plot(T, nCa17/nCa, label = 'nCa17', color = '#f7b6d2')
plt.plot(TEvolx, nCa17x/nCax, label = 'nCa17x', color = '#f7b6d2', linestyle = ':')
plt.plot(T, nCa18/nCa, label = 'nCa18', color = '#c7c7c7')
plt.plot(TEvolx, nCa18x/nCax, label = 'nCa18x', color = '#c7c7c7', linestyle = ':')
plt.plot(T, nCa19/nCa, label = 'nCa19', color = '#dbdb8d')
plt.plot(TEvolx, nCa19x/nCax, label = 'nCa19x', color = '#dbdb8d', linestyle = ':')
plt.plot(T, nCa20/nCa, label = 'nCa20', color = '#9edae5')
plt.plot(TEvolx, nCa20x/nCax, label = 'nCa20x', color = '#9edae5', linestyle = ':')
plt.yscale('log')
plt.xscale('log')
plt.ylim(2e-3, 1.2)
plt.xlim(1e4, 1e7)
plt.savefig("nCa_vs_time.png", dpi = 300)
plt.close()

plt.plot(T, nFe0/nFe, label = 'nFe0', color = '#1f77b4')
plt.plot(TEvolx, nFe0x/nFex, label = 'nFe0x', color = '#1f77b4', linestyle = ':')
plt.plot(T, nFe1/nFe, label = 'nFe1', color = '#ff7f0e')
plt.plot(TEvolx, nFe1x/nFex, label = 'nFe1x', color = '#ff7f0e', linestyle = ':')
plt.plot(T, nFe2/nFe, label = 'nFe2', color = '#2ca02c')
plt.plot(TEvolx, nFe2x/nFex, label = 'nFe2x', color = '#2ca02c', linestyle = ':')
plt.plot(T, nFe3/nFe, label = 'nFe3', color = '#d62728')
plt.plot(TEvolx, nFe3x/nFex, label = 'nFe3x', color = '#d62728', linestyle = ':')
plt.plot(T, nFe4/nFe, label = 'nFe4', color = '#9467bd')
plt.plot(TEvolx, nFe4x/nFex, label = 'nFe4x', color = '#9467bd', linestyle = ':')
plt.plot(T, nFe5/nFe, label = 'nFe5', color = '#8c564b')
plt.plot(TEvolx, nFe5x/nFex, label = 'nFe5x', color = '#8c564b', linestyle = ':')
plt.plot(T, nFe6/nFe, label = 'nFe6', color = '#e377c2')
plt.plot(TEvolx, nFe6x/nFex, label = 'nFe6x', color = '#e377c2', linestyle = ':')
plt.plot(T, nFe7/nFe, label = 'nFe7', color = '#7f7f7f')
plt.plot(TEvolx, nFe7x/nFex, label = 'nFe7x', color = '#7f7f7f', linestyle = ':')
plt.plot(T, nFe8/nFe, label = 'nFe8', color = '#bcbd22')
plt.plot(TEvolx, nFe8x/nFex, label = 'nFe8x', color = '#bcbd22', linestyle = ':')
plt.plot(T, nFe9/nFe, label = 'nFe9', color = '#17becf')
plt.plot(TEvolx, nFe9x/nFex, label = 'nFe9x', color = '#17becf', linestyle = ':')
plt.plot(T, nFe10/nFe, label = 'nFe10', color = '#1f77b4')
plt.plot(TEvolx, nFe10x/nFex, label = 'nFe10x', color = '#1f77b4', linestyle = ':')
plt.plot(T, nFe11/nFe, label = 'nFe11', color = '#aec7e8')
plt.plot(TEvolx, nFe11x/nFex, label = 'nFe11x', color = '#aec7e8', linestyle = ':')
plt.plot(T, nFe12/nFe, label = 'nFe12', color = '#ffbb78')
plt.plot(TEvolx, nFe12x/nFex, label = 'nFe12x', color = '#ffbb78', linestyle = ':')
plt.plot(T, nFe13/nFe, label = 'nFe13', color = '#98df8a')
plt.plot(TEvolx, nFe13x/nFex, label = 'nFe13x', color = '#98df8a', linestyle = ':')
plt.plot(T, nFe14/nFe, label = 'nFe14', color = '#ff9896')
plt.plot(TEvolx, nFe14x/nFex, label = 'nFe14x', color = '#ff9896', linestyle = ':')
plt.plot(T, nFe15/nFe, label = 'nFe15', color = '#c5b0d5')
plt.plot(TEvolx, nFe15x/nFex, label = 'nFe15x', color = '#c5b0d5', linestyle = ':')
plt.plot(T, nFe16/nFe, label = 'nFe16', color = '#c49c94')
plt.plot(TEvolx, nFe16x/nFex, label = 'nFe16x', color = '#c49c94', linestyle = ':')
plt.plot(T, nFe17/nFe, label = 'nFe17', color = '#f7b6d2')
plt.plot(TEvolx, nFe17x/nFex, label = 'nFe17x', color = '#f7b6d2', linestyle = ':')
plt.plot(T, nFe18/nFe, label = 'nFe18', color = '#c7c7c7')
plt.plot(TEvolx, nFe18x/nFex, label = 'nFe18x', color = '#c7c7c7', linestyle = ':')
plt.plot(T, nFe19/nFe, label = 'nFe19', color = '#dbdb8d')
plt.plot(TEvolx, nFe19x/nFex, label = 'nFe19x', color = '#dbdb8d', linestyle = ':')
plt.plot(T, nFe20/nFe, label = 'nFe20', color = '#9edae5')
plt.plot(TEvolx, nFe20x/nFex, label = 'nFe20x', color = '#9edae5', linestyle = ':')
plt.plot(T, nFe21/nFe, label = 'nFe21', color = '#393b79')
plt.plot(TEvolx, nFe21x/nFex, label = 'nFe21x', color = '#393b79', linestyle = ':')
plt.plot(T, nFe22/nFe, label = 'nFe22', color = '#5254a3')
plt.plot(TEvolx, nFe22x/nFex, label = 'nFe22x', color = '#5254a3', linestyle = ':')
plt.plot(T, nFe23/nFe, label = 'nFe23', color = '#6b6ecf')
plt.plot(TEvolx, nFe23x/nFex, label = 'nFe23x', color = '#6b6ecf', linestyle = ':')
plt.plot(T, nFe24/nFe, label = 'nFe24', color = '#9c9ede')
plt.plot(TEvolx, nFe24x/nFex, label = 'nFe24x', color = '#9c9ede', linestyle = ':')
plt.plot(T, nFe25/nFe, label = 'nFe25', color = '#637939')
plt.plot(TEvolx, nFe25x/nFex, label = 'nFe25x', color = '#637939', linestyle = ':')
plt.plot(T, nFe26/nFe, label = 'nFe26', color = '#8ca252')
plt.plot(TEvolx, nFe26x/nFex, label = 'nFe26x', color = '#8ca252', linestyle = ':')
plt.yscale('log')
plt.xscale('log')
plt.ylim(2e-3, 1.2)
plt.xlim(1e4, 1e7)
plt.savefig("nFe_vs_time.png", dpi = 300)
plt.close()

print('Done !!!')

print('Elapsed time = ', time.time() - TA)

