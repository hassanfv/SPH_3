
import numpy as np
import matplotlib.pyplot as plt
import pickle
import glob


pc_to_cm = 3.086e18

#===== closestNdx
def closestNdx(arr, val):
  return np.argmin(abs(arr - val))

dirX = '/home/pc/Desktop/N_body_2024/SPH_3/cooling29June2024/Table_preparation_3/On_EC2_Oct12/pklFiles01/pklFiles/'

filez = glob.glob(dirX + '*.pkl')

#----------- Preparing the grid -------
rkpcG = np.arange(0.01, 1.02, 0.1)# it is in kpc
LshG = np.arange(0.0, 2.51, 0.25) # it is in log10 of pc so 0.0 mean 1pc or 3.086e18 cm ---> We take Lsh in the range 1.0 pc up to ~300 pc.
nHG = np.arange(-4.0, 4.01, 0.1)
#--------------------------------------

nH = 3.2
rkpc = 0.21
Lsh = 0.75

ndx_nH = closestNdx(nHG, nH)
ndx_rkpc = closestNdx(rkpcG, rkpc)
ndx_Lsh = closestNdx(LshG, Lsh)

nam = dirX + f'./nH_{nHG[ndx_nH]:.1f}_rkpc_{rkpcG[ndx_rkpc]:.2f}_Lsh_{np.log10(10**LshG[ndx_Lsh] * pc_to_cm):.3f}.pkl'
print(nam)

with open(nam, 'rb') as f:
  data = pickle.load(f)

TEvol = data['TempEvol']
t_Arr_in_yrs = data['t_in_sec'] / 3600. / 24. / 365.25


#----- Plot Section --------
plt.scatter(t_Arr_in_yrs, np.log10(TEvol), s = 1)

plt.show()
#--------------------------








