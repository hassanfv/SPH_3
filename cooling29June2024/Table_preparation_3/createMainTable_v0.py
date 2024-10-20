
import numpy as np
import matplotlib.pyplot as plt
import pickle
import glob

#===== closestNdx
def closestNdx(arr, val):
  return np.argmin(abs(arr - val))


#===== closestNdx
def boundNdx(arr, val):
  nx = np.argmin(abs(arr - val))
  if val >= arr[nx]:
    nxLeft = nx
    nxRight= nx + 1
  else:
    nxLeft = nx - 1
    nxRight = nx
  return nxLeft, nxRight


pc_to_cm = 3.086e18

#----------- Preparing the grid -------
rkpcG = np.arange(0.01, 1.02, 0.1)# it is in kpc
LshG = np.arange(0.0, 2.51, 0.25) # it is in log10 of pc so 0.0 mean 1pc or 3.086e18 cm ---> We take Lsh in the range 1.0 pc up to ~300 pc.
nHG = np.arange(-4.0, 4.01, 0.1)
TG = np.arange(2.0, 10.31, 0.1)
#--------------------------------------

#-------- Creating mu grid ------------
n_steps = int((1.25 - 0.59) / (np.mean([0.04, 0.01])))
print('len n_steps = ', n_steps)
muSteps = np.linspace(0.04, 0.01, n_steps)
y0 = 0.59
muG = np.zeros_like(muSteps)
for i, tmp in enumerate(muSteps):
  muG[i] = y0
  y0 += tmp
#--------------------------------------

print()
print('muG = ', muG)
print('len(muG) = ', len(muG))



dirX = '/home/pc/Desktop/N_body_2024/SPH_3/cooling29June2024/Table_preparation_3/pklFromEC2/'
filez = glob.glob(dirX + '*.pkl')
print(len(filez))
print()

j = 7

nam = filez[j]

print(nam)

with open(nam, 'rb') as f:
  data = pickle.load(f)
# 'TempEvol', 'AbundEvol', 'nH', 'rkpc', 'Lsh', 'Species_id', 'Species_name', 't_in_sec', 'nH_p', 'rkpc_p', 'Lsh_p', 'mu'
keyz = list(data.keys())

TEvol = np.log10(data['TempEvol'])
AbEvol = data['AbundEvol']
Species_id = data['Species_id']
Species_name = data['Species_name']
nH_p = float(data['nH_p'])
rkpc_p = float(data['rkpc_p'])
Lsh_p = float(data['Lsh_p'])
t_Arr_in_yrs = data['t_in_sec'] / 3600. / 24. / 365.25

print()
print(f'nH_p, rkpc_p, Lsh_p = {nH_p}, {rkpc_p}, {Lsh_p}')

print()
print(keyz)
print()

if 'mu' in keyz:
  muEvol = data['mu']
else:
  muEvol = 0.6 + np.zeros_like(TEvol)

print('muEvol = ', muEvol)
print()

#----- Finding the indices for nH, rkpc, Lsh, and mu to be used in mainTable -------
ndx_nH = closestNdx(nHG, nH_p)
ndx_rkpc = closestNdx(rkpcG, rkpc_p)
ndx_Lsh = closestNdx(LshG, Lsh_p)

#ndx_Lsh = closestNdx(muG, mu_p) # mu_p changes over time but other _p parameters are fixed over time!!!!
#-----------------------------------------------------------------------------------

# Before assigning mu_p I need to find the location of T_i for interpolation!!!!!! 

print('TEvol = ', TEvol)

T_p = 5.4

nx = closestNdx(TEvol, T_p)

T1 = TEvol[nx-1]
T2 = TEvol[nx]
T3 = TEvol[nx+1]
Tarr = [T1, T2, T3]

t1 = t_Arr_in_yrs[nx-1]
t2 = t_Arr_in_yrs[nx]
t3 = t_Arr_in_yrs[nx+1]
tarr = [t1, t2, t3]
tZeroPoint = min(tarr)

t1 -= tZeroPoint
t2 -= tZeroPoint
t3 -= tZeroPoint

delta_T = np.abs(T_p - T2)
print(f'\ndelta_T = {delta_T}\n')

if delta_T < 0.01:
  #doInterpCase1(t1, t2, t3, T1, T2, T3)
  
  tarr = [t1, t2, t3]
  Tarr = [T1, T2, T3]
  
  t100s = t2 + np.arange(0, 101)
  
  T_interp = np.interp(t100s, tarr, Tarr)
  
  #print(f'\nT_interp = {T_interp}\n')
  

#if delta_T >= 0.01:




#new_t = np.linspace(min(tarr), max(tarr), 1000)
#T_interp = np.interp(new_t, tarr, Tarr)

#nx2 = closestNdx(T_interp, T_p)

print('\n---------------')
print(f'First closest T to {T_p} is {T2}\n')
#print(f'Second closest T to {T_p} is {T_interp[nx2]}\n')
print('T1, T2, T3 = ', T1, T2, T3)
print()
print('t after zeroPoint: ', t1, t2, t3)
print()
print(tZeroPoint)
print('\n')
#print('T_interp = ', T_interp)

s()


ndx_T_L, ndx_T_R = boundNdx(TEvol, T_p)

print()
print(TEvol[ndx_T_L], T_p, TEvol[ndx_T_R])

s()

nH = 2.5
rkpc = 0.31
#Lsh = 0.75
Lsh = np.log10(10**19.489 / 3.086e18)

ndx_nH = closestNdx(nHG, nH)
ndx_rkpc = closestNdx(rkpcG, rkpc)
ndx_Lsh = closestNdx(LshG, Lsh)

nam = dirX + f'./nH_{nHG[ndx_nH]:.1f}_rkpc_{rkpcG[ndx_rkpc]:.2f}_Lsh_{np.log10(10**LshG[ndx_Lsh] * pc_to_cm):.3f}.pkl'
print(nam)

with open(nam, 'rb') as f:
  data = pickle.load(f)

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





