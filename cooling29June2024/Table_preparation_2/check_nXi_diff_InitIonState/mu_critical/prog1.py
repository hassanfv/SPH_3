import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import h5py
import os


df = pd.read_csv('data_species.csv')
print(df)
    
AtomicMass = df['A']

gamma = 5./3.
kB = 1.3807e-16
mH = 1.6726e-24
pc_to_cm = 3.086e18


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

  user_input = {
    "output_file": "   " + OutFile + '\n',
    "distance_to_AGN_kpc": f'{rkpc_i:.2f}\n',
    "max_shield_length": " " + f'{(10**Lsh_i * pc_to_cm):.3E}\n',
    "log_nH_min": 4*" " + f'{log_nH_i:.1f}\n',
    "log_nH_max": 4*" " + f'{log_nH_i:.1f}\n'}

  # Update the parameters in the file content
  updated_content = update_parameters(original_content, user_input)

  # Write the updated content to a new file
  updated_file_name = f'nH_{log_nH_i:.1f}_rkpc_{rkpc_i:.2f}_Lsh_{np.log10(10**Lsh_i * pc_to_cm):.3f}.param'
  with open(updated_file_name, 'w') as file:
      file.writelines(updated_content)

  #---- Executing CHIMES ----
  command = f"python3 chimes-driver.py {updated_file_name}"
  #command = f"mpirun -np 96 -x LD_LIBRARY_PATH=/path/to/install/dir/lib:$LD_LIBRARY_PATH python3 chimes-driver.py {updated_file_name}"
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

  N_nH = n_densities = f['TableBins/N_Densities'][()]
  N_T = n_temperatures = f['TableBins/N_Temperatures'][()]
  nHarr = f['TableBins/Densities'][:]
  Tarr = f['TableBins/Temperatures'][:]
  
  return TempEvol, AbundEvol, nHarr, Tarr





#----------- Preparing the grid -------
rkpcG = np.arange(0.01, 1.02, 0.1)# it is in kpc
LshG = np.arange(0.0, 2.51, 0.25) # it is in log10 of pc so 0.0 mean 1pc or 3.086e18 cm ---> We take Lsh in the range 1.0 pc up to ~300 pc.
nHG = np.arange(1.0, 4.01, 0.1)
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

ResX = np.zeros((N, 2))

for j in range(N):

  nH_t, rkpc_t, Lsh_t = inList[j]

  OutHDF5FileName, OutFile_pkl, updated_file_name = mainFunc(nH_t, rkpc_t, Lsh_t)
  
  print()
  print('OutHDF5FileName, OutFile_pkl = ', OutHDF5FileName, OutFile_pkl)
  
  TempEvol, AbundEvol, nHarr, Tarr = getModelOutput(OutHDF5FileName)
  
  TempEvol = TempEvol[0, 0, 0, :]
    
  print()
  print('TempEvol.shape = ', TempEvol.shape)
  print()
  
  AbundEvol = AbundEvol[0, 0, 0, :, :]
  
  print()
  print('AbundEvol.shape = ', AbundEvol.shape)
  print()
  
  #------ Finding the turning point -----
  Tmp1 = TempEvol[:-1]
  Tmp2 = TempEvol[1:]

  nxx = np.where(((Tmp1 - Tmp2) < 0.0) & (Tmp2 < 1.1e4))[0] # Assuming turning point always occurs below 1e4 K. 

  print('nxx = ', nxx)

  if nxx.size > 0:
    ndx = nxx[0]
    T_crit = TempEvol[ndx]
    AbX = AbundEvol[:, ndx]
    mu_crit = getMu(AbX)
    ResX[j, 0] = T_crit
    ResX[j, 1] = mu_crit
    print('ndx = ', ndx)
    print('T_crit = ', T_crit)
    print('mu_crit = ', mu_crit)
  else:
    print('nxx is empty !!!!!')
  
  os.remove(OutHDF5FileName)
  os.remove(updated_file_name)
  
  print(ResX)

  #s()




