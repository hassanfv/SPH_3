
import numpy as np
import pandas as pd
import pickle
import os
import time
import h5py
import subprocess


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


#----------- Preparing the grid -------
rkpcG = np.arange(0.01, 1.02, 0.2)# it is in kpc
LshG = np.arange(0.0, 2.51, 0.5) # it is in log10 ---> We take Lsh in the range 1.0 pc up to ~300 pc.
#rkpcG = np.arange(0.01, 1.02, 0.2)# it is in kpc
#LshG = np.arange(0.0, 2.51, 0.5) # it is in log10 ---> We take Lsh in the range 1.0 pc up to ~300 pc.
#--------------------------------------


#===== mainFunc
def mainFunc(i, rkpc_i, Lsh_i):

  OutFile = f'./grid_{i:03}_rkpc_{rkpc_i:.2f}_Lsh_{np.log10(10**Lsh_i * pc_to_cm):.3f}.hdf5'

  user_input = {
    "output_file": "   " + OutFile + '\n',
    "distance_to_AGN_kpc": f'{rkpc_i:.3f}\n',
    "max_shield_length": " " + f'{(10**Lsh_i * pc_to_cm):.3E}\n'}

  # Update the parameters in the file content
  updated_content = update_parameters(original_content, user_input)

  # Write the updated content to a new file
  updated_file_name = f'grid_{i:03}_rkpc_{rkpc_i:.2f}_Lsh_{np.log10(10**Lsh_i * pc_to_cm):.3f}.param'
  with open(updated_file_name, 'w') as file:
      file.writelines(updated_content)

  #---- Executing CHIMES ----
  #command = f"mpirun -np 1 python3 chimes-driver.py {updated_file_name}"
  command = f"mpirun -np 5 -x LD_LIBRARY_PATH=/path/to/install/dir/lib:$LD_LIBRARY_PATH python3 chimes-driver.py {updated_file_name}"
  print(command)
  os.system(command)
  #----
  
  os.remove(updated_file_name)

  return OutFile


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



pc_to_cm = 3.086e18

SpName = ['HI', 'OI', 'CI', 'CII', 'CIV', 'SiII', 'SiIV', 'NV', 'FeII', 'MgI', 'MgII', 'SII']
#                HI OI CI CII CIV SiII SiIV  NV  FeII  MgI  MgII  SII
SelectSpecies = [1, 23, 7, 8, 10, 58,  60,   19, 111,  44,  45,   73]

# Read the initial param file.
with open('Refgrid_noneq_evolution_AGN.param', 'r') as file:
    original_content = file.readlines()


TA = time.time()

N_Species = len(SelectSpecies)
N_rkpc = len(rkpcG)
N_Lsh = len(LshG)
N_T = 37 # calculate by looking at grid_***.param !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
N_nH = 41 # calculate by looking at grid_***.param !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
N_time = 101 # calculate by looking at grid_***.param !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

TEvol = np.zeros((N_rkpc, N_Lsh, N_nH, N_T, N_time))
AbEvol = np.zeros((N_rkpc, N_Lsh, N_nH, N_T, N_Species, N_time))

print('N_rkpc * N_Lsh = ', N_rkpc * N_Lsh)

counter = 0
for i in range(N_rkpc):
  for j in range(N_Lsh):
  
    OutHDF5FileName = mainFunc(counter, rkpcG[i], LshG[j])
    counter += 1
    
    TempEvol, AbundEvol, nHarr, Tarr = getModelOutput(OutHDF5FileName)
    
    print()
    print('TempEvol.shape = ', TempEvol.shape)
    print()
    
    TB = time.time()
    for k in range(N_nH):
      for l in range(N_T):
        TEvol[i, j, k, l, :] = np.log10(TempEvol[l, k, 0, :]) # TEvol[rkpc, Lsh, nH, T, time]
        for b in range(N_Species):
          tmp = np.log10(1e-30+AbundEvol[l, k, 0, SelectSpecies[b], :])
          nt = np.where(tmp < -30.0)[0]
          if len(nt) > 0:
            tmp[nt] = -30.0
          AbEvol[i, j, k, l, b, :] = tmp
    
    os.remove(OutHDF5FileName)


TEvolx = {'TEvol': TEvol, 'nH': nHarr, 'Temp': Tarr, 'rkpc': rkpcG, 'Lsh': LshG, 't': np.arange(N_time)}
with open('TEvol.pkl', 'wb') as f:
  pickle.dump(TEvolx, f)

AbEvolx = {'AbEvol': AbEvol, 'nH': nHarr, 'Temp': Tarr, 'rkpc': rkpcG, 'Lsh': LshG, 'Species_id': SelectSpecies, 'Species_name': SpName, 't': np.arange(N_time)}
with open('AbEvol.pkl', 'wb') as f:
  pickle.dump(AbEvolx, f)


print(f'Elapsed timeX = {(time.time() - TA):.1f} sec')





