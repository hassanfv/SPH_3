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


if False:
  print('----> reactants <----')
  print(reactants[:, :])
  print()
  print('----> products <----')
  print(products[:, :])
  print()


#print(np.sort(reactants[:, 0]))


nt = np.where(reactants[:, 0] == 1)[0]   # Change the index here to see the reactions and products of that element
print('******** reactants ****************')
print('nt = ', nt)
print(reactants[nt, :])
print('************************')
print()

print('******** products ****************')
print(products[nt, :])
print('************************')
print()

N = len(nt)

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
  
  N = 31 - len(strx)
  strx = strx + N * ' ' + f'   (ndx = {nt[j]})'
  

  #label = ' ----> charge transfer'
  #if 'γ' in strx:
  #  label = ' ----> radiative recombination'
  #if 'e + e' in strx:
  #  label = ' ----> collisional ionization'
  #strz = strx + label
  
  print(strx)







