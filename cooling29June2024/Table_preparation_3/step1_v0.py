import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import h5py
import os
from scipy.interpolate import interp1d


#=========================================
df = pd.read_csv('nH_vs_time.csv')
# ['nH', 'n_iter', 'dt']
print(df.columns)

nH = df['nH'].values
n_iter = df['n_iter'].values
dt = df['dt'].values

t = n_iter * dt
t_yrs = t / 3600. / 24. / 365.25

#--- Interpolation ------
nHG = np.linspace(-4.0, 4.0, 81)  # Generates 81 points from -4 to 4 inclusively
f_linear = interp1d(nH, np.log10(t_yrs), bounds_error=False, fill_value="extrapolate")
log_tG = f_linear(nHG)

t_G = 10**log_tG * 3600. * 24. * 365.25 # converting back to seconds!
dt_G = t_G / 30000.0

print(dt_G)
#=========================================


df = pd.read_csv('xdata_species.csv')
print(df)
    
AtomicMass = df['A']

gamma = 5./3.
kB = 1.3807e-16
mH = 1.6726e-24
pc_to_cm = 3.086e18

SpName = ['HI', 'OI', 'CI', 'CII', 'CIV', 'SiII', 'SiIV', 'NV', 'FeII', 'MgI', 'MgII', 'SII']
#                HI OI CI CII CIV SiII SiIV  NV  FeII  MgI  MgII  SII
SelectSpecies = [1, 23, 7, 8, 10, 58,  60,   19, 111,  44,  45,   73]


# Read the initial param file.
with open('Refgrid_noneq_evolution_AGN.param', 'r') as file:
    original_content = file.readlines()


#----- getMu
def getMu(Ab):

  s = 0.0
  p = 0.0
  for j in range(157):
    s += Ab[j] * AtomicMass[j]
    p += Ab[j] # Note that ne is also included in the sum!!

  mu = s / p
  
  return mu


#===== update_parameters
def update_parameters(content, updates):
    
    for i, line in enumerate(content):
        # Split the line into components assuming space or tab delimitation
        parts = line.split()
        # Check if the line contains a parameter that needs to be updated
        if len(parts) > 1 and parts[0] in updates:
            # Update the parameter value based on the user input
            parts[1] = str(updates[parts[0]])
            # Reconstruct the line with the updated value
            content[i] = ' '.join(parts)
    return content



#===== mainFunc
def mainFunc(log_nH_i, rkpc_i, Lsh_i):

  OutFile = f'./nH_{log_nH_i:.1f}_rkpc_{rkpc_i:.2f}_Lsh_{np.log10(10**Lsh_i * pc_to_cm):.3f}.hdf5'
  OutFile_pkl = f'./nH_{log_nH_i:.1f}_rkpc_{rkpc_i:.2f}_Lsh_{np.log10(10**Lsh_i * pc_to_cm):.3f}.pkl'
  
  ndx = (np.abs(nHG - log_nH_i)).argmin()
  dt_i = dt_G[ndx]

  user_input = {
    "output_file": "   " + OutFile + '\n',
    "distance_to_AGN_kpc": f'{rkpc_i:.2f}\n',
    "max_shield_length": " " + f'{(10**Lsh_i * pc_to_cm):.3E}\n',
    "log_nH_min": 4*" " + f'{log_nH_i:.1f}\n',
    "log_nH_max": 4*" " + f'{log_nH_i:.1f}\n',
    "hydro_timestep": 11*" " + f'{dt_i:.2E}\n'}

  # Update the parameters in the file content
  updated_content = update_parameters(original_content, user_input)

  # Write the updated content to a new file
  updated_file_name = f'nH_{log_nH_i:.1f}_rkpc_{rkpc_i:.2f}_Lsh_{np.log10(10**Lsh_i * pc_to_cm):.3f}.param'
  with open(updated_file_name, 'w') as file:
      file.writelines(updated_content)

  #---- Executing CHIMES ----
  command = f"python3 chimes-driver.py {updated_file_name}"
  print(command)
  os.system(command)
  #----
  
  #os.remove(updated_file_name)

  return OutFile, OutFile_pkl, updated_file_name # updated_file_name will be used to delete param file at the end of each loop!


#===== getModelOutput
def getModelOutput(OutHDF5FileName):

  f = h5py.File(OutHDF5FileName, 'r')
  TempEvol = f['TemperatureEvolution'][:] # (T, nH, Z, t)
  AbundEvol = f['AbundanceEvolution'][:]     # (T, nH, Z, Elm, t)
  t_Arr_in_sec = f['TimeArray_seconds'][:]

  N_nH = n_densities = f['TableBins/N_Densities'][()]
  N_T = n_temperatures = f['TableBins/N_Temperatures'][()]
  nHarr = f['TableBins/Densities'][:]
  #Tarr = f['TableBins/Temperatures'][:]
  
  return TempEvol, AbundEvol, t_Arr_in_sec





#----------- Preparing the grid -------
rkpcG = np.arange(0.01, 1.02, 0.1)# it is in kpc
LshG = np.arange(0.0, 2.51, 0.25) # it is in log10 of pc so 0.0 mean 1pc or 3.086e18 cm ---> We take Lsh in the range 1.0 pc up to ~300 pc.
nHG = np.arange(-4.0, 4.01, 0.1)
#--------------------------------------

print(len(rkpcG), len(LshG), len(nHG), len(rkpcG)*len(LshG)*len(nHG))

N_nH = len(nHG)
N_rkpc = len(rkpcG)
N_Lsh = len(LshG)

inList = [] #np.zeros(N_nH * N_rkpc * N_Lsh)

for i in range(N_nH):
  for j in range(N_rkpc):
    for k in range(N_Lsh):
      inList.append([float(nHG[i]), float(rkpcG[j]), float(LshG[k])])


N = len(inList)

#ResX = np.zeros((N, 2))

for j in range(N):

  nH_t, rkpc_t, Lsh_t = inList[j]

  OutHDF5FileName, OutFile_pkl, updated_file_name = mainFunc(nH_t, rkpc_t, Lsh_t)
  
  print()
  print('OutHDF5FileName, OutFile_pkl = ', OutHDF5FileName, OutFile_pkl)
  
  TempEvol, AbundEvol, t_Arr_in_sec = getModelOutput(OutHDF5FileName)
  
  TempEvol = TempEvol[0, 0, 0, :]
    
  print()
  print('TempEvol.shape = ', TempEvol.shape)
  print()
  
  
  AbundEvol = AbundEvol[0, 0, 0, SelectSpecies, :]
    
  print()
  print('AbundEvol.shape = ', AbundEvol.shape)
  print()
  
  dictx = {'TempEvol': TempEvol, 'AbundEvol': AbundEvol,
           'nH': nHG, 'rkpc': rkpcG,
           'Lsh': LshG, 'Species_id': SelectSpecies,
           'Species_name': SpName, 't_in_sec': t_Arr_in_sec}
  
  with open(OutFile_pkl, 'wb') as f:
    pickle.dump(dictx, f)
  
  #os.remove(OutHDF5FileName)
  #os.remove(updated_file_name)
  
  #print(ResX)

  #s()




