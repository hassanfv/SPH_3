import h5py
import numpy as np
import matplotlib.pyplot as plt


# Reaction: (Hp + e ---> H0 + γ)  ::: photo-recombination of Hp and e.
def k2(T):

  k2_val = 8.4e-11 / T**0.5 / (T/1e3)**0.2 / (1.0 + (T/1e6)**0.7)

  return k2_val

# Reaction: (Hep + e ---> He0 + γ)  ::: photo-recombination of Hep and e.
def k5(T):

  k5_val = 1.5e-10 / T**0.6353
  
  return k5_val
#--------------

#Reaction: Hep di-electric recombination (Hep + e ---> He0 + γ)
def k7(T):

  k7_val = 1.9e-3 / T**1.5 * np.exp(-470000./T) * (1.0 + 0.3 * np.exp(-94000./T))
  
  return k7_val
#--------------

#-------> cooling rates from Cen - 1992 <-----------
# Cooling rate due to H collisional ionization: (H0 + e ---> Hp + 2e)  ::: in erg.s^-1.cm^3  ::: will be multiplied by nH0 and ne later in the code.
def g1(T):

  g1_val = 1.27e-21 * T**0.5 * np.exp(-157809.1/T) / (1.0 + (T/1e5)**0.5)
  
  return g1_val


# Cooling rate due to Hp photo-recombination: (Hp + e ---> H0 + γ)  ::: in erg.s^-1.cm^3  ::: will be multiplied by nHp and ne later in the code.
def g2(T):

  #g2_val = 8.70e-27 * T**0.5 / (T/1e3)**0.2 / (1.0 + (T/1e6)**0.7)
  g2_val = kB * T * k2(T) # See LaMothe and Ferland - 2001, page 1, the lines below equation 1!

  return g2_val


# Cooling due to H0 collisional excitation via collision with electrons ::: in erg.s^-1.cm^3 ::: will be multiplied by nH0 and ne later in the code.
def g3(T):

  g3_val = 7.50e-19 * np.exp(-118348/T) / (1.0 + (T/1e5)**0.5)
  
  return g3_val


# Cooling due to Hep collisional excitation ::: will be multiplied by nHep and ne later in the code
def g5(T):

  g5_val = 5.54e-17 / T**0.397 * np.exp(-473638./T) / (1. + (T/1e5)**0.5)
  
  return g5_val


# Cooling due to He0 collisional ioniozation ::: will be multiplied by nHe0 and ne later in the code
def g6(T):
  
  g6_val = 9.38e-22 * T**0.5 * np.exp(-285335.4/T) / (1. + (T/1e5)**0.5)
  
  return g6_val


# Cooling via Hep collisional ionization ::: will be multiplied by nHep and ne later in the code
def g7(T):

  g7_val = 4.95e-22 * T**0.5 * np.exp(-631515./T) / (1. + (T/1e5)**0.5)
  
  return g7_val


# Cooling via Hep recombination with electron ::: will be multiplied by nHep and ne later in the code
def g8(T):

  g8_val = 1.55e-26 * T**0.3647
  
  return g8_val


# Cooling via Hepp recombination with electron ::: will be multiplied by nHepp and ne later in the code
def g9(T):
  
  g9_val = 3.48e-26 * T**0.5 / (T/1e3)**0.2 / (1. + (T/1e6)**0.7)
  
  return g9_val


# Cooling via di-electronic recombination of Hep with electrons ::: will be multiplied by nHep and ne later in the code
def g10(T):
  
  g10_val = 1.24e-13 / T**1.5 * np.exp(-470000./T) * (1. + 0.3 * np.exp(-94000./T))

  return g10_val



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

rateA_HII = rateCaseA[0, :]
rateA_HeII = rateCaseA[1, :]

rateB_HII = rateCaseB[0, :]
rateB_HeII = rateCaseB[1, :]


#-------> Cent - 1992 <-------
r_HII = np.log10(k2(10**Temp))
r_HeII = np.log10(k5(10**Temp + k7(10**Temp)))
#-----------------------------


#-----------------------------------------------------------
#------------> Cooling rates from CHIMES table <------------
#-----------------------------------------------------------
with h5py.File('chimes_main_data.hdf5', 'r') as file:

  coolants = file['cooling/coolants'][:]
  cooling_rates = file['cooling/rates'][:]

print()
print('coolants.shape: ', coolants.shape)
print(coolants)
print()
print('cooling_rates.shape: ', cooling_rates.shape)
print()
print(cooling_rates)
print()
# We plot cooling by HI, HII, HeI, HeII, HeIII or indices: 0, 1, 2, 3, 4

G_HI    = 10**cooling_rates[0, :]
G_HII   = 10**cooling_rates[1, :]
G_HeI   = 10**cooling_rates[2, :]
G_HeII  = 10**cooling_rates[3, :]
G_HeIII = 10**cooling_rates[4, :]

G_H = np.log10(G_HI + G_HII)
G_He= np.log10(G_HeI + G_HeII + G_HeIII)



kB = 1.3807e-16 # erg/K
gH = np.log10(g1(10**Temp) + g2(10**Temp) + g3(10**Temp))
gHe= np.log10(g5(10**Temp) + g6(10**Temp) + g7(10**Temp) + g8(10**Temp) + g9(10**Temp) + g10(10**Temp))


plt.figure(figsize = (14, 6))

plt.subplot(1, 2, 1)
plt.scatter(Temp, rateA_HII, color = 'k', s = 3, label = 'HII + e ---> HI + 2e')
plt.scatter(Temp, rateA_HeII, color = 'b',s = 3, label = 'HeII + e ---> HeI + 2e')

plt.plot(Temp, rateB_HII, color = 'k', label = 'Case B - HII + e ---> HI + 2e', linestyle = '--')
plt.plot(Temp, rateB_HeII, color = 'b', label = 'Case B - HeII + e ---> HeI + 2e', linestyle = '--')

#plt.scatter(Temp, r_HII, color = 'orange', s = 3, label = 'RR_HII from Cen - 1992')
#plt.scatter(Temp, r_HeII, color = 'lime', s = 3, label = 'RR_HeII from Cen - 1992')
plt.legend()


plt.subplot(1, 2, 2)
plt.scatter(Temp, G_H, color = 'k', s = 3, label = 'cooling via H')
plt.scatter(Temp, G_He, color = 'b',s = 3, label = 'cooling via He')

plt.scatter(Temp, gH, color = 'orange', s = 3, label = 'Cen-1992 - cooling via H')
plt.scatter(Temp, gHe, color = 'lime', s = 3, label = 'Cen-1992 - cooling via He')

plt.legend()


plt.show()



