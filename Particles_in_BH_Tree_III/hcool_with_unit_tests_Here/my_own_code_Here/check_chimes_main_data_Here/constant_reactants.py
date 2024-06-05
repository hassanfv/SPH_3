import h5py
import numpy as np


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






# Open the HDF5 file
with h5py.File('chimes_main_data.hdf5', 'r') as file:
  
  N_reactions = file['constant/N_reactions'][:]
  print('N_reactions = ', N_reactions)


  reactants = file['constant/reactants'][:]
  products = file['constant/products'][:]


print("Shape of 'constant/reactants':", reactants.shape)
print("Data type of 'constant/reactants':", reactants.dtype)
print(np.sort(reactants[:, 0]))
print()
print('******* All reactants ******')
print((reactants))
print()
print()
print('+++++++ All products +++++++')
print(products)
print()


#------ Printing the reactants and products with their elemental names for better visualization -------
a = reactants[:, 0]
b = reactants[:, 1]

x = products[:, 0]
y = products[:, 1]
z = products[:, 2]

for i in range(len(a)):

  if y[i] != -1 and z[i] != -1:
    print(f'{elmList[a[i]]} + {elmList[b[i]]} ---> {elmList[x[i]]} + {elmList[y[i]]} + {elmList[z[i]]}')
  
  if y[i] != -1 and z[i] == -1:
    print(f'{elmList[a[i]]} + {elmList[b[i]]} ---> {elmList[x[i]]} + {elmList[y[i]]}')
  
  if y[i] == -1 and z[i] == -1:
    print(f'{elmList[a[i]]} + {elmList[b[i]]} ---> {elmList[x[i]]}')


#------------------------------------------------------------------------------------------

s()


i = 2

print(reactants[i, :])

N = np.sum(reactants[i, :] != -1)

for j in range(N):
  print(f'Reactant_{i} = {elmList[reactants[i, j]]}')

print()
nt = np.where(reactants[:, 0] == 3)[0]   # Change the index here to see the reactions and products of that element
print('******** reactants ****************')
print('nt = ', nt)
print(reactants[nt, :])
print('************************')
print()

print('******** products ****************')
print(products[nt, :])
print('************************')






