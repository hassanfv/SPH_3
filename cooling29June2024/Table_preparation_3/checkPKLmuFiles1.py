
import numpy as np
import matplotlib.pyplot as plt
import pickle
import glob


pc_to_cm = 3.086e18

#===== closestNdx
def closestNdx(arr, val):
  return np.argmin(abs(arr - val))

#dirX = '/home/pc/Desktop/N_body_2024/SPH_3/cooling29June2024/Table_preparation_3/pklMuFiles/'
#filez = glob.glob(dirX + '*.pkl')

#----------- Preparing the grid -------
rkpcG = np.arange(0.01, 1.02, 0.1)# it is in kpc
LshG = np.arange(0.0, 2.51, 0.25) # it is in log10 of pc so 0.0 mean 1pc or 3.086e18 cm ---> We take Lsh in the range 1.0 pc up to ~300 pc.
nHG = np.arange(-4.0, 4.01, 0.1)
#--------------------------------------

nH = 2.5
rkpc = 0.31
#Lsh = 0.75
Lsh = np.log10(10**19.489 / 3.086e18)

ndx_nH = closestNdx(nHG, nH)
ndx_rkpc = closestNdx(rkpcG, rkpc)
ndx_Lsh = closestNdx(LshG, Lsh)

#nam = dirX + f'./nH_{nHG[ndx_nH]:.1f}_rkpc_{rkpcG[ndx_rkpc]:.2f}_Lsh_{np.log10(10**LshG[ndx_Lsh] * pc_to_cm):.3f}.pkl'
#print(nam)

dirX = '/home/pc/Desktop/N_body_2024/SPH_3/cooling29June2024/Table_preparation_3/pklFromEC2/'
filez = glob.glob(dirX + '*.pkl')
j = 5771
nam = filez[j]

with open(nam, 'rb') as f:
  data = pickle.load(f)
nH_p = float(data['nH_p'])
rkpc_p = float(data['rkpc_p'])
Lsh_p = float(data['Lsh_p'])

print()
print(f'nH_p, rkpc_p, Lsh_p = {nH_p}, {rkpc_p}, {Lsh_p}\n')

print()
print(data.keys())
print()

TEvol = data['TempEvol']
t_Arr_in_yrs = data['t_in_sec'] / 3600. / 24. / 365.25

mu = data['mu']

#----- Plot Section --------
#plt.scatter(t_Arr_in_yrs, np.log10(TEvol), s = 1)
plt.scatter(np.log10(TEvol), mu, s = 1)

plt.show()
#--------------------------





