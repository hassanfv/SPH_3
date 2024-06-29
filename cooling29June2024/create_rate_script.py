import h5py
import numpy as np


elmList = [
            "e", "HI", "HII", "Hm", "HeI", "HeII", "HeIII", "CI", "CII", "CIII",
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


Mol = ["H2", "H2p", "H3p", "OH", "H2O", "C2", "O2", "HCOp", "CH",
       "CH2", "CH3p", "CO", "CHp", "CH2p", "OHp", "H2Op", "H3Op", "COp", "HOCp", "O2p"]


res1_list = []
res2_list = []

# Open the HDF5 file
with h5py.File('chimes_main_data.hdf5', 'r') as file:

  N_reactions = file['T_dependent/N_reactions'][:]
  print('N_reactions = ', N_reactions)
  print()
  
  reactants = file['T_dependent/reactants'][:]
  print("reactants.shape':", reactants.shape)
  products = file['T_dependent/products'][:]
  print("products.shape':", reactants.shape)
  print()


oneTimerX = 0

for i in range(len(elmList)):

  nt = np.where(reactants[:, 0] == i)[0]

  N = len(nt)

  a = b = c = x = y = z = ''

  for j in range(N):

    reac = reactants[nt[j], :]

    #===== The reactants section ======
    a = elmList[reac[0]]
    b = elmList[reac[1]]  
    if reac[2] != -1:
      c = elmList[reac[2]]
    
    str1 = f'{a} + {b}'
    if reac[2] != -1:
      str1 = str1 + f' + {c}'

    str1 = str1 + f' ----> '

    #===== The products section ======
    prod = products[nt[j], :]

    x = elmList[prod[0]]
    if prod[1] != -1:
      y = elmList[prod[1]]
    if prod[2] != -1:
      z = elmList[prod[2]]

    str2 = f'{x}'
    if prod[1] != -1:
      str2 = str2 + f' + {y}'
    if prod[2] != -1:
      str2 = str2 + f' + {z}'
    
    if str2 == f'{x}':
      str2 = str2 + f' + γ'
    
    strx = str1 + str2
    
    spaces = ' ' * (31 - len(strx))
    
    N = 31 - len(strx)
    strx = strx + N * ' ' + f'   (ndx = {nt[j]})'
    
    label = ''
    
    if 'γ' in strx:
      label = ' ----> radiative recombination'
    if 'e + e' in strx:
      label = ' ----> collisional ionization'
    
    if (a in Mol) or (b in Mol) or (c in Mol) or (x in Mol) or (y in Mol) or (z in Mol):
      label = '----> MOLECULES involved !!!'
    
    strz = strx + label
    
    # ['H', 'He', 'C', 'N', 'O', 'Ne', 'Mg', 'Si', 'S', 'Ca', 'Fe']
    
    if label != '----> MOLECULES involved !!!' and a[0] == 'F': #!!!!!!!!!!!! CONFLICT with (C and Ca) AND (N and Ne) ???????!!!!!! Note a[0] is only one 
                                                                # charactor so to get Ne and Ca copy them from the output for N and C ! clever !!!
                                                                # For Mg only type 'M', for Si type S, for Fe type F !clever!
      
      res1 = f'R_{a}_to_{x}_via_{b}_ = rates[{nt[j]}, :]'
      res2 = f'R_{a}_to_{x}_via_{b} = interp1d(Temp, R_{a}_to_{x}_via_{b}_, kind="linear", fill_value="extrapolate")'
      
      if (a[0] == 'H') and (oneTimerX == 0):
        tmp = 'R_HII_to_HI_via_e_caseA_ = ratesAB[0, :] # H CaseA'
        res1_list.append(tmp)
        tmp = 'R_HeII_to_HeI_via_e_caseA_ = ratesAB[1, :] # He CaseA'
        res1_list.append(tmp)
        tmp = 'R_HII_to_HI_via_e_caseA = interp1d(Temp, R_HII_to_HI_via_e_caseA_, kind="linear", fill_value="extrapolate") # H CaseA'
        res2_list.append(tmp)
        tmp = 'R_HeII_to_HeI_via_e_caseA = interp1d(Temp, R_HeII_to_HeI_via_e_caseA_, kind="linear", fill_value="extrapolate") # He CaseA'
        res2_list.append(tmp)
        oneTimerX = 1
      
      res1_list.append(res1)
      res2_list.append(res2)
      
      print(res1)
      print(res2)

# Writing to text files
with open('res1.py', 'w') as file1:
    for item in res1_list:
        file1.write(item + "\n")

with open('res2.py', 'w') as file2:
    for item in res2_list:
        file2.write(item + "\n")


print()
print('-----------------------------------------------------------------')
print('-------------- Results saved to res1.py and res2.py -------------')
print('-----------------------------------------------------------------')
print()


