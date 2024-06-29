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

#----- getAtmNum
def getAtmNum(iD):
  iDlist = np.array(['C', 'N', 'O', 'Ne', 'Mg', 'Si', 'S', 'Ca', 'Fe'])
  AtmNumlist = [6, 7, 8, 10, 12, 14, 16, 20, 26]
  n = np.where(iD == iDlist)[0][0]
 
  return AtmNumlist[n]


elm = 'Fe'   # !!!!!!!!!!!!!!!!!!!!!! Change this for other elements and re-run the code !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
AtmNum = getAtmNum(elm)

spec_list = [elm+roman_num[i] for i in range(AtmNum+1)] # current species list

print(spec_list)
print()


#-----------------------------------------------------------
#------------> Cooling rates from CHIMES table <------------
#-----------------------------------------------------------
with h5py.File('chimes_main_data.hdf5', 'r') as file:

  coolants = file['cooling/coolants'][:]
  cooling_rates = file['cooling/rates'][:]

print()
print('coolants.shape: ', coolants.shape)
print()
print('cooling_rates.shape: ', cooling_rates.shape)
print()

print('++++++++ coolants ++++++++')
print(coolants)
print()
print()

for j in range(len(coolants)):
  print(f'{coolants[j]} --> {elmList[coolants[j]]}')
print()


res1 = []
res2 = []

res1_list = []
res2_list = []


for j in range(len(coolants)):
  
  if elmList[coolants[j]] in spec_list:
    #print(j, coolants[j], elmList[coolants[j]])
    res1 = f'g{elmList[coolants[j]]}_ = cooling_rates[{j}, :]'
    res2 = f'g{elmList[coolants[j]]} = interp1d(Temp, g{elmList[coolants[j]]}_, kind="linear", fill_value="extrapolate")'
      
    res1_list.append(res1)
    res2_list.append(res2)
    
    print(res1)
    print(res2)

# Writing to text files
with open('lmb1.py', 'w') as file1:
    for item in res1_list:
        file1.write(item + "\n")

with open('lmb2.py', 'w') as file2:
    for item in res2_list:
        file2.write(item + "\n")

print()
print('-----------------------------------------------------------------')
print('-------------- Results saved to lmb1.py and lmb2.py -------------')
print('-----------------------------------------------------------------')
print()






