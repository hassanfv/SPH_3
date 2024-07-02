import h5py
import numpy as np
import matplotlib.pyplot as plt


elmList = [
            "elec", "HI", "HII", "Hm", "HeI", "HeII", "HeIII", "CI", "CII", "CIII",
            "CIV", "CV", "CVI", "CVII", "Cm", "NI", "NII", "NIII", "NIV", "NV",
            "NVI", "NVII", "NVIII", "OI", "OII", "OIII", "OIV", "OV", "OVI", "OVII",
            "OVIII", "OIX", "Om", "NeI", "NeII", "NeIII", "NeIV", "NeV", "NeVI",
            "NeVII", "NeVIII", "NeIX", "NeX", "NeXI", "MgI", "MgII", "MgIII", "MgIV",
            "MgV", "MgVI", "MgVII", "MgVIII", "MgIX", "MgX", "MgXI", "MgXII", "MgXIII",
            "SiI", "SiII", "SiIII", "SiIV", "SiV", "SiVI", "SiVII", "SiVIII", "SiIX",
            "SiX", "SiXI", "SiXII", "SiXIII", "SiXIV", "SiXV", "SI", "SII", "SIII",
            "SIV", "SV", "SVI", "SVII", "SVIII", "SIX", "SX", "SXI", "SXII", "SXIII",
            "SXIV", "SXV", "SXVI", "SXVII", "CaI", "CaII", "CaIII", "CaIV", "CaV",
            "CaVI", "CaVII", "CaVIII", "CaIX", "CaX", "CaXI", "CaXII", "CaXIII", "CaXIV",
            "CaXV", "CaXVI", "CaXVII", "CaXVIII", "CaXIX", "CaXX", "CaXXI", "FeI",
            "FeII", "FeIII", "FeIV", "FeV", "FeVI", "FeVII", "FeVIII", "FeIX", "FeX",
            "FeXI", "FeXII", "FeXIII", "FeXIV", "FeXV", "FeXVI", "FeXVII", "FeXVIII",
            "FeXIX", "FeXX", "FeXXI", "FeXXII", "FeXXIII", "FeXXIV", "FeXXV", "FeXXVI",
            "FeXXVII", "H2", "H2p", "H3p", "OH", "H2O", "C2", "O2", "HCOp", "CH",
            "CH2", "CH3p", "CO", "CHp", "CH2p", "OHp", "H2Op", "H3Op", "COp", "HOCp", "O2p"
          ]


roman_num = [
                  "I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X",
                  "XI", "XII", "XIII", "XIV", "XV", "XVI", "XVII", "XVIII", "XIX", "XX",
                  "XXI", "XXII", "XXIII", "XXIV", "XXV", "XXVI", "XXVII"
                 ]

#neutrals = ['CI', 'NI', 'OI', 'NeI', 'MgI', 'SiI', 'SI', 'CaI', 'FeI']

excludeList = ['CI', 'CII', 'NII', 'OI', 'SiII', 'FeII', 'Cm', 'Om'] # Except Cm and Om, others are treated as 2d and 4d !! double check to make sure!!

grain_recomb_list = ['CII', 'OII', 'SiII', 'FeII', 'MgII', 'SII', 'CaII', 'CaIII'] # HII and HeII are not included as they are hard-coded!

#----- getAtmNum
def getAtmNum(iD):
  iDlist = np.array(['C', 'N', 'O', 'Ne', 'Mg', 'Si', 'S', 'Ca', 'Fe'])
  AtmNumlist = [6, 7, 8, 10, 12, 14, 16, 20, 26]
  n = np.where(iD == iDlist)[0][0]
 
  return AtmNumlist[n]


elm = 'C'   
AtmNum = getAtmNum(elm)
spec_list = [elm+roman_num[i] for i in range(AtmNum+1)] # current species list
spec_iState = [i for i in range(AtmNum+1)] # current species ionization states
spec_list += ['Cm']
spec_iState += [-1]

if True:
  #---------------------------------------------------------------------------------------------------
  elm = 'N'   
  AtmNum = getAtmNum(elm)
  spec_list2 = [elm+roman_num[i] for i in range(AtmNum+1)] # current species list
  spec_iState2 = [i for i in range(AtmNum+1)] # current species ionization states
  spec_list += spec_list2
  spec_iState += spec_iState2
  #----------------------------------------------------------------------------------------------------

if True:
  #---------------------------------------------------------------------------------------------------
  elm = 'O'   
  AtmNum = getAtmNum(elm)
  spec_list2 = [elm+roman_num[i] for i in range(AtmNum+1)] # current species list
  spec_iState2 = [i for i in range(AtmNum+1)] # current species ionization states
  spec_list += spec_list2
  spec_iState += spec_iState2
  spec_list += ['Om']
  spec_iState += [-1]
  #---------------------------------------------------------------------------------------------------
  
if True:
  #---------------------------------------------------------------------------------------------------
  elm = 'Ne'   
  AtmNum = getAtmNum(elm)
  spec_list2 = [elm+roman_num[i] for i in range(AtmNum+1)] # current species list
  spec_iState2 = [i for i in range(AtmNum+1)] # current species ionization states
  spec_list += spec_list2
  spec_iState += spec_iState2
  
if True:
  #---------------------------------------------------------------------------------------------------
  elm = 'Mg'   
  AtmNum = getAtmNum(elm)
  spec_list2 = [elm+roman_num[i] for i in range(AtmNum+1)] # current species list
  spec_iState2 = [i for i in range(AtmNum+1)] # current species ionization states
  spec_list += spec_list2
  spec_iState += spec_iState2

if True:
  #---------------------------------------------------------------------------------------------------
  elm = 'Si'   
  AtmNum = getAtmNum(elm)
  spec_list2 = [elm+roman_num[i] for i in range(AtmNum+1)] # current species list
  spec_iState2 = [i for i in range(AtmNum+1)] # current species ionization states
  spec_list += spec_list2
  spec_iState += spec_iState2

if True:
  #---------------------------------------------------------------------------------------------------
  elm = 'S'   
  AtmNum = getAtmNum(elm)
  spec_list2 = [elm+roman_num[i] for i in range(AtmNum+1)] # current species list
  spec_iState2 = [i for i in range(AtmNum+1)] # current species ionization states
  spec_list += spec_list2
  spec_iState += spec_iState2

if True:
  #---------------------------------------------------------------------------------------------------
  elm = 'Ca'   
  AtmNum = getAtmNum(elm)
  spec_list2 = [elm+roman_num[i] for i in range(AtmNum+1)] # current species list
  spec_iState2 = [i for i in range(AtmNum+1)] # current species ionization states
  spec_list += spec_list2
  spec_iState += spec_iState2

if True:
  #---------------------------------------------------------------------------------------------------
  elm = 'Fe'   
  AtmNum = getAtmNum(elm)
  spec_list2 = [elm+roman_num[i] for i in range(AtmNum+1)] # current species list
  spec_iState2 = [i for i in range(AtmNum+1)] # current species ionization states
  spec_list += spec_list2
  spec_iState += spec_iState2

print(spec_list)
print(spec_iState)
print()



strx = '#----- Lambda\n'
strx += f'def Lambda(T, nHI, nHII, nHm, nHeI, nHeII, nHeIII, '

k = 6 # to take into account the already existance of nHp, nHe, nHepp!
for x in spec_list:
  strx += f'n{x}, '
  if not(k % 10) and x != spec_list[-1]:
    strx += f'\n{11 * " "}'
  k += 1

strx = strx[:-2]

#if 'CI' in spec_list: # 'CI' is only used just to know that Carbon is activated... we could use any other Carbon ions name here!
#  strx += ', nCm'

#if 'OI' in spec_list:
#  strx += ', nOm'

strx += '):\n\n'


strx +=f'  Tx = np.log10(T)\n\n'

str_ne = f'  ne = (\n{8 * " "} 1 * nHII - nHm + (nHeII + 2.0 * nHeIII)'

kk = 3
for k, x in zip(spec_iState, spec_list):
  if k !=0:
    str_ne += f' + {k} * n{x}'
    if not (kk % 6):
      str_ne += f'\n{6 * " "}'
  kk += 1

str_ne += f'\n{7 * " "})\n\n'
strx += str_ne


strx += '  #----- # Glover & Jappsen - 2007 -----\n'
strx += f'  z = 0.0 # current time redshift!\n'
strx += f'  TCMB_0 = 2.7255\n'
strx += f'  TCMB = TCMB_0 * (1.0 + z)\n'
strx += f'  LCompton = 1.017e-37 * TCMB**4 * (T - TCMB) * ne\n'
strx += f'  #--------------------------------------\n\n'

# NOTE: grain cooling rate is the same for all species, HII, HeII, CII, OII, SiII, FeII, MgII, SII, CaII, CaIII--> (See cooling/grain_recombination)
strx += f'  grain_recomb_cooling_rate = grain_cool_rate(T, ne, nHI, nHII, G0, A_v, dust_ratio, Temp, Psi)\n\n'
 


strx += f'  Lamb = (\n'
strx += f'  {7 * " "} 10**g1(Tx) * ne * nHI  # H0\n'
strx += f'  {5 * " "} + 10**g2(Tx) * ne * nHII # Hp\n'
strx += f'  {5 * " "} + 10**grain_recomb_cooling_rate * nHII * ne # grain_recombination cooling!\n' 
strx += f'  {5 * " "} + 10**g3(Tx) * nHeI * ne # He0\n' 
strx += f'  {5 * " "} + 10**g4(Tx) * nHeII * ne # Hep\n'
strx += f'  {5 * " "} + 10**grain_recomb_cooling_rate * nHeII * ne # grain_recombination cooling!\n' 
strx += f'  {5 * " "} + 10**g5(Tx) * nHeIII * ne# Hep\n'

if 'CI' in spec_list: # 'CI' is only used just to know that Carbon is activated... we could use any other Carbon ions name here!
  strx += f'  {5 * " "} + 10**cooling_rate_4d("CI", T, nHI, ne, nHII, Temp_4d, HIDensity_4d, elecDensity_4d, HIIDensity_4d) * nCI * ne # cooling via CI\n'
  strx += f'  {5 * " "} + 10**cooling_rate_2d("CII", T, ne, Temp_2d, elecDensity_2d) * nCII * ne # cooling via CII\n'
if 'NII' in spec_list:
  strx += f'  {5 * " "} + 10**cooling_rate_2d("NII", T, ne, Temp_2d, elecDensity_2d) * nNII * ne # cooling via NII\n'
if 'OI' in spec_list:
  strx += f'  {5 * " "} + 10**cooling_rate_4d("OI", T, nHI, ne, nHII, Temp_4d, HIDensity_4d, elecDensity_4d, HIIDensity_4d) * nOI * ne # cooling via OI\n'
if 'SiII' in spec_list:
  strx += f'  {5 * " "} + 10**cooling_rate_2d("SiII", T, ne, Temp_2d, elecDensity_2d) * nSiII * ne # cooling via SiII\n'
if 'FeII' in spec_list:
  strx += f'  {5 * " "} + 10**cooling_rate_2d("FeII", T, ne, Temp_2d, elecDensity_2d) * nFeII * ne # cooling via FeII\n'

for x in spec_list:
  print('x = ', x)
  if x not in excludeList:
    strx += f'  {5 * " "} + 10**g{x}(Tx) * n{x} * ne\n'
  if x in grain_recomb_list:
    strx += f'  {5 * " "} + 10**grain_recomb_cooling_rate * n{x} * ne # grain_recombination cooling!\n' 

strx += f'  {5 * " "} + LCompton)\n\n'
#strx += f'  {5 * " "} - (Lcr_H0 + Lcr_He0 + Lcr_Hep + Lcr_C0 + Lcr_C1 + Lcr_C2 + Lcr_C3 + Lcr_C4 + Lcr_C5)) # (-) multiplied by (-) become (+).\n\n'

strx += f'  return Lamb'

print(strx)

with open('output.py', 'w') as file:
    file.write(strx)

print('\n\n')
print('------------------------------------------------')
print()
print('        File saved to output.py file       ')
print()
print('------------------------------------------------')
print()


