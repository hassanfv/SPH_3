
import numpy as np
import pandas as pd
import pickle
import os
import time
import h5py
import subprocess
import glob


# Now we remove hdf5 files only after the inner loop is fully completed (check line "hdf5_files = glob.glob("*.hdf5")")! (06 September - 2024)
# I added: "Checking if hdf5 file is created and visible to the system if not wait for it to happen" (06 September - 2024)

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
rkpcG = np.arange(0.01, 0.52, 0.1)# it is in kpc
LshG = np.arange(0.0, 2.51, 0.25) # it is in log10 ---> We take Lsh in the range 1.0 pc up to ~300 pc.
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
  command = f"mpirun -np 96 -x LD_LIBRARY_PATH=/path/to/install/dir/lib:$LD_LIBRARY_PATH python3 chimes-driver.py {updated_file_name}"
  print(command)
  os.system(command)
  #----
  
  #os.remove(updated_file_name)

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
N_T = 73 # calculate by looking at grid_***.param !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
N_nH = 81 # calculate by looking at grid_***.param !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
N_time = 101 # calculate by looking at grid_***.param !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

TEvol = np.zeros((N_rkpc, N_Lsh, N_nH, N_T, N_time))
AbEvol = np.zeros((N_rkpc, N_Lsh, N_nH, N_T, N_Species, N_time))

print('N_rkpc * N_Lsh = ', N_rkpc * N_Lsh)

counter = 0
for i in range(N_rkpc):
  for j in range(N_Lsh):
  
    TAA = time.time()
    
    debuger = 0
  
    while debuger == 0:
  
      try:
        OutHDF5FileName = mainFunc(counter, rkpcG[i], LshG[j])
        
        #------ Checking if hdf5 file is created and visible to the system if not wait for it to happen ------
        file_exists = False
        timeout = 10  # seconds
        start_time = time.time()

        while not file_exists and (time.time() - start_time) < timeout:
          file_exists = os.path.exists(OutHDF5FileName)
          if not file_exists:
            time.sleep(1)  # Sleep for a second and then check again

        if not file_exists:
          print("HASSAN: Failed to find hdf5 file after waiting!!!!!!!!!!!!!!!!!!!!!!!!!!")
        else:
          print("HASSAN: hdf5 File found, proceeding!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        #------------------------------------------------------------------------------------------------------
        
        TempEvol, AbundEvol, nHarr, Tarr = getModelOutput(OutHDF5FileName)
        debuger = 1 # if hdf5 file exists then debuger will be 1 and it will exit the while loop !!
        counter += 1
      except:
        print('\n\n\n')
        print('EXCEPTION OCCURED !!!!!!!!!!!!!!! PROBLEM --> hdf5 file can not be created by chimes !!!!!!!!!!!!!!')
        print('\n\n\n')
      
    
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

    print(f'Elapse time in one main loop = {time.time() - TAA} seconds')
    print()

  # Loop through the files and remove them
  '''
  hdf5_files = glob.glob("*.hdf5")
  for filex in hdf5_files:
    try:
      os.remove(filex)
      print(f"Removed {filex}")
    except Exception as e:
      print(f"Error removing {filex}: {e}")
  '''
    


TEvolx = {'TEvol': TEvol, 'nH': nHarr, 'Temp': Tarr, 'rkpc': rkpcG, 'Lsh': LshG, 't': np.arange(N_time)}
with open('TEvol.pkl', 'wb') as f:
  pickle.dump(TEvolx, f)

AbEvolx = {'AbEvol': AbEvol, 'nH': nHarr, 'Temp': Tarr, 'rkpc': rkpcG, 'Lsh': LshG, 'Species_id': SelectSpecies, 'Species_name': SpName, 't': np.arange(N_time)}
with open('AbEvol.pkl', 'wb') as f:
  pickle.dump(AbEvolx, f)


print(f'Elapsed timeX = {(time.time() - TA):.1f} sec')





