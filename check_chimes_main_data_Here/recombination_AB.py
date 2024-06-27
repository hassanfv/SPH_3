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



# Open the HDF5 file
with h5py.File('chimes_main_data.hdf5', 'r') as file:

  N_reactions = file['recombination_AB/N_reactions'][:]
  print('N_reactions = ', N_reactions)

  reactants = file['recombination_AB/reactants'][:]
  products = file['recombination_AB/products'][:]
  
  rateCaseA = file['recombination_AB/rates_caseA'][:]
  print('rateCaseA.shape', rateCaseA.shape)
  
  rateCaseB = file['recombination_AB/rates_caseB'][:]
  print('rateCaseB.shape', rateCaseB.shape)
  
  Temp = file['TableBins/Temperatures'][:]
  print('Temp.shape', Temp.shape)
  print()
  
  element_incl = file['recombination_AB/element_incl'][:]
  print('element_incl:')
  print(element_incl)
  print()
  
  molecular_flag = file['recombination_AB/molecular_flag'][:]
  print('molecular_flag:', molecular_flag)
  print()


print("Shape of 'recombination_AB/reactants':", reactants.shape)
print("Data type of 'recombination_AB/reactants':", reactants.dtype)
print()
print('---- reactants ----')
print('2 is HII, 5 is HeII and 0 is electron.')
print('So here we see collision of HII and HeII with electron (radiative recombination)')
print(reactants[:, :])
print('-------------------')
print()
print('---- products ----')
print('1 is HI    4 is HeI')
print('So the recombination of HII and HeII with electron results in HI and HeI')
print(products[:])
print('-------------------')




