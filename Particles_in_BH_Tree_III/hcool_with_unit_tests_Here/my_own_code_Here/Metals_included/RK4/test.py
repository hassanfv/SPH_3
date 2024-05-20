
import numpy as np
import pandas as pd


#----------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
#------------------------------------- Cooling Rates (erg.s^-1.cm^3) --------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
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


# Cooling due to free-free emission ::: will be multiplied by (nHp+nHep+nHepp)*ne later in the code.
def g4(T):

  gff = 1.1 + 0.34 * np.exp(-(5.5 - np.log10(T))**2/3.0)
  g4_val = 1.42e-27 * gff * T**0.5
  
  return g4_val


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


df = pd.read_csv('CData.csv')
Tref = df['Temp'].values


Tgrid = np.logspace(3, 9, 1000)


g1x = np.interp(Tref, x_known, y_known)







