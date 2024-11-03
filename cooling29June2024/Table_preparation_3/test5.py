
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle
import glob


#----------- Preparing the grid -------
rkpcG = np.arange(0.01, 1.02, 0.1)# it is in kpc
LshG = np.arange(0.0, 2.51, 0.25) # it is in log10 of pc so 0.0 mean 1pc or 3.086e18 cm ---> We take Lsh in the range 1.0 pc up to ~300 pc.
nHG = np.arange(-4.0, 4.01, 0.1)
TG = np.arange(2.0, 10.31, 0.1)

#-------- Creating mu grid ------------
#n_steps = int((1.25 - 0.59) / (np.mean([0.04, 0.01])))
#print('len n_steps = ', n_steps)
#muSteps = np.linspace(0.04, 0.01, n_steps)
#y0 = 0.59
#muG = np.zeros_like(muSteps)
#for i, tmp in enumerate(muSteps):
#  muG[i] = y0
#  y0 += tmp

#muG = np.arange(0.6, 1.26, 0.05)

muG = np.array([0.6, 0.7, 0.8, 0.9, 0.95, 1.0] + [float(_) for _ in np.arange(1.01, 1.24, 0.01)])

N_mu = len(muG)
#--------------------------------------

print(muG)
print('N_mu = ', N_mu)



N_nH = len(nHG)
N_rkpc = len(rkpcG)
N_Lsh = len(LshG)
N_T = len(TG)
N_mu = len(muG)

#------ TEST FILE ----------------------
dirX = '/home/pc/Desktop/N_body_2024/SPH_3/cooling29June2024/Table_preparation_3/pklFromEC2/'
filez = glob.glob(dirX + '*.pkl')
print(len(filez))
print()
j = 5771  #196#11
nam = filez[j]
with open(nam, 'rb') as f:
  data = pickle.load(f)
nH_i = nH_p = float(data['nH_p'])
rkpc_i = rkpc_p = float(data['rkpc_p'])
Lsh_i = Lsh_p = float(data['Lsh_p'])
muEvol = data['mu']
print(f'nH_p = {nH_p},  rkpc_p = {rkpc_p},  Lsh_p = {Lsh_p}\n')


TEvol = np.log10(data['TempEvol'])
t_Arr_in_yrs = data['t_in_sec'] / 3600. / 24. / 365.25
#---------------------------------------

with open('mainTable.pkl', 'rb') as f:
  data = pickle.load(f)

with open('mainTableMu.pkl', 'rb') as f:
  dataMu = pickle.load(f)
# mainTable = np.zeros((N_nH, N_rkpc, N_Lsh, N_mu, N_T, 101))

#--- T_i = 7.6 is our starting point
T_i = 7.0 # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
mu_i = 0.6 # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#Lsh_i = 0.0
# Step 1: Find the index fo, nHi, rkpc_i, Lsh_i, mu_i, and T_i in nHG, rkpcG, LshG, muG, and TG !
# NOTE: We also need the evolved mu, i.e. the mu after the dt timestep evolution !

t = 0
res = []
res.append([t, T_i, mu_i])

debugRes = []


for j in range(1000):

  nx_nH_i = round((nH_i - nHG[0]) / 0.1)
  nx_rkpc_i = round((rkpc_i - rkpcG[0]) / 0.1)
  nx_Lsh_i = round((Lsh_i - LshG[0]) / 0.25)
  nx_mu_i = np.argmin(abs(muG - mu_i))
  nx_T_i = round((T_i - TG[0]) / 0.1)

  TEvolved = data[nx_nH_i, nx_rkpc_i, nx_Lsh_i, nx_mu_i, nx_T_i][:15]   #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  MuEvolved = dataMu[nx_nH_i, nx_rkpc_i, nx_Lsh_i, nx_mu_i, nx_T_i][:15]#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ', TG[nx_T_i])
  
  if TG[nx_T_i] <= 5.300009:
    print('Above duck-tapping Section -------')
    print(f'len(TEvolved) = {len(TEvolved)},    len(MuEvolved) = {len(MuEvolved)}')
    print('TEvolved = ', TEvolved)
    print()
    print('MuEvolved = ', MuEvolved)
    print()
    print('---------------------------------\n')
  
  #------------ duck-tapping ------
  k = 0
  TEvolved = data[nx_nH_i, nx_rkpc_i, nx_Lsh_i, nx_mu_i+k, nx_T_i][:15]
  while (TEvolved[0] == 0) & (nx_mu_i+k < N_mu-1):
    
    print('mu in LOOP = ', muG[nx_mu_i+k])
    TEvolved = data[nx_nH_i, nx_rkpc_i, nx_Lsh_i, nx_mu_i+k, nx_T_i][:15]   #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    MuEvolved = dataMu[nx_nH_i, nx_rkpc_i, nx_Lsh_i, nx_mu_i+k, nx_T_i][:15]#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    k += 1
  
  if k > 0 or True:
    print(f'ATTEN.: k = {k} so it is greater than 0 !!!!')
    print('j = ', j)
    
    #for kk in range(-1, N_mu-1):
    #  if (nx_mu_i+kk) < 26:
    #    TEvolved1 = data[nx_nH_i, nx_rkpc_i, nx_Lsh_i, nx_mu_i+kk, nx_T_i][:]
    #    print(f'k = {float(kk)}:  {int(nx_mu_i),   float(muG[nx_mu_i+kk]), float(TG[nx_T_i]), float(TEvolved1[0])}')

  delta_T = TEvolved[-1] - TEvolved[0]
  
  delta_Mu= MuEvolved[-1] - MuEvolved[0]

  T_i_old = T_i    # for debugging!
  mu_i_old = mu_i  # for debugging!

  T_i += delta_T
  mu_i += delta_Mu
  
  t += len(TEvolved)
  
  res.append([t, T_i, mu_i])
  
  print(T_i_old, TG[nx_T_i], TEvolved[0], TEvolved[-1], delta_T, T_i, mu_i_old, muG[nx_mu_i], k, MuEvolved[0], MuEvolved[-1], delta_Mu, mu_i)
  
  '''
  if TG[nx_T_i] <= 5.300009:
    print(f'k = ', k)
    print(f'len(TEvolved) = {len(TEvolved)},    len(MuEvolved) = {len(MuEvolved)}')
    print('TEvolved = ', TEvolved)
    print()
    print('MuEvolved = ', MuEvolved)
    print()
    s()
  '''
  
  debugRes.append([T_i_old, TG[nx_T_i], TEvolved[0], TEvolved[-1], delta_T, T_i, mu_i_old, muG[nx_mu_i], k, MuEvolved[0], MuEvolved[-1], delta_Mu, mu_i])
  

res = np.array(res)
t = res[:, 0]
T = res[:, 1]
mu = res[:, 2]

dfBug = pd.DataFrame(debugRes)
dfBug.columns = ['T_i_old', 'T_closest', 'TStart', 'TEnd', 'd_T', 'T_new', 'mu_i_old', 'mu_closest', 'k', 'muStart', 'muEnd', 'd_mu', 'mu_new']
dfBug.to_csv('deBug.csv', index = False)


#--- Test plot ---
# Create a figure with two subplots side by side
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 5))

# First scatter plot
ax1.scatter(t_Arr_in_yrs, TEvol, color='k', s=20)
ax1.scatter(t+952757, T, color='r', s=5)
ax1.set_title('TEvol vs Time')
ax1.set_xlabel('t_Arr_in_yrs')
ax1.set_ylabel('TEvol')

# Second scatter plot
ax2.scatter(t_Arr_in_yrs, muEvol, color='k', s=20)
ax2.scatter(t+952757, mu, color='r', s=5)
ax2.set_title('muEvol vs Time')
ax2.set_xlabel('t_Arr_in_yrs')
ax2.set_ylabel('muEvol')

# Third scatter plot
ax3.scatter(TEvol, muEvol, color='k', s=20)
ax3.scatter(T, mu, color='r', s=5)
ax3.set_title('muEvol vs T')
ax3.set_xlabel('T')
ax3.set_ylabel('mu')

# Display the plots
plt.tight_layout()
plt.show()
#-----------------



#print(tEvolved)
#print()
#print(MuEvolved)

# Now update T_i and mu_i and go to the next round. BTW, keep track of time, current T_i and mu_i for later plotting !!!




