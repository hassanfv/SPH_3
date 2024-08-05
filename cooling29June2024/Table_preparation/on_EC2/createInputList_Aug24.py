
import numpy as np
import pickle


rkpc = np.arange(0.02, 0.73, 0.1)
logLsh = np.arange(12, 23.5, 1.0) # Lsh = NHtot / nH... with 16 < NHtot < 23 and -4 < nH < 4 we get 12 < Lsh < 23.5. Note 23.5 is log(3.086e23)==log(max_Lsh)

print(logLsh)

N_rkpc = len(rkpc)
N_Lsh = len(logLsh)

print('rkpc = ', rkpc)
print()

print(f'N_rkpc = {N_rkpc}, N_Lsh = {N_Lsh}')

N_tot = N_rkpc * N_Lsh

print(f'Number of CHIMES models to be generated = {N_tot}')

res = []

for k in range(N_rkpc):
  for l in range(N_Lsh):
    res.append([rkpc[k], logLsh[l]]) # Note that we take logLsh to the power of 10 to make it cm !

res = np.array(res)
print()
print(res.shape)

with open('inputLists_Aug24.pkl', 'wb') as f:
  pickle.dump(res, f)



