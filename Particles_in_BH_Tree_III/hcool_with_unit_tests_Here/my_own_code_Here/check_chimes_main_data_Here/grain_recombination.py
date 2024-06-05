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


# Only for these species we have grain_recombination reaction!!!
nt = [2, 5, 8, 24, 58, 111, 45, 73, 90, 91]
#   2      5     8     24    58     111    45     73    90      91
#['HII' 'HeII' 'CII' 'OII' 'SiII' 'FeII' 'MgII' 'SII' 'CaII' 'CaIII']

print(np.array(elmList)[nt])
print()


# Open the HDF5 file
with h5py.File('chimes_main_data.hdf5', 'r') as file:

  N_reactions = file['grain_recombination/N_reactions'][:]
  print('N_reactions = ', N_reactions)
  print()
  
  reactants = file['grain_recombination/reactants'][:]
  print("reactants.shape':", reactants.shape)
  products = file['grain_recombination/products'][:]
  print("products.shape':", reactants.shape)
  print()
  
  grain_recomn_rates = file['grain_recombination/rates'][:]
  print('grain_recomn_rates.shape = ', grain_recomn_rates.shape)
  
  Temp = file['TableBins/Temperatures'][:]
  print('Temp.shape', Temp.shape)
  print()
  
  Psi = file['TableBins/Psi'][:]
  print('Psi.shape', Temp.shape)
  print()
  
  grain_cooling_rates = file['cooling/grain_recombination'][:] 
  print('grain_cooling_rates.shape = ', grain_cooling_rates.shape)
  print()


if True:
  print('----> reactants <----')
  print(reactants[:, :])
  print()
  print('----> products <----')
  print(products)
  print()


aa = grain_recomn_rates[0, :]
print(aa.shape)



