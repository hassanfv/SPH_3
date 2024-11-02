
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle


#----------- Preparing the grid -------
rkpcG = np.arange(0.01, 1.02, 0.1)# it is in kpc
LshG = np.arange(0.0, 2.51, 0.25) # it is in log10 of pc so 0.0 mean 1pc or 3.086e18 cm ---> We take Lsh in the range 1.0 pc up to ~300 pc.
nHG = np.arange(-4.0, 4.01, 0.1)
TG = np.arange(2.0, 10.31, 0.1)


with open('mainTable.pkl', 'rb') as f:
  data = pickle.load(f)

# mainTable = np.zeros((N_nH, N_rkpc, N_Lsh, N_mu, N_T, 101))

i, j, k, l = [75, 6, 0, 24] # nH, rkpc, Lsh, mu

tmp = data[i, j, k, 0, 20] # Note to change mu when playing with T.

print(tmp)
print()
print(tmp.shape)
print(TG[70])
print()
print(nHG[i])




