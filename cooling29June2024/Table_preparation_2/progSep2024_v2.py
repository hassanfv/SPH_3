
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
rkpcG = np.arange(0.01, 1.02, 0.1)# it is in kpc
LshG = np.arange(0.0, 2.51, 0.25) # it is in log10 ---> We take Lsh in the range 1.0 pc up to ~300 pc.

lst = []
for rkpc in rkpcG:
  for Lsh in LshG:
    lst.append([rkpc, Lsh])
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
  command = f"mpirun -np 5 python3 chimes-driver.py {updated_file_name}"
  print(command)
  os.system(command)
  #----

  return OutFile


pc_to_cm = 3.086e18

# Read the initial param file.
with open('Refgrid_noneq_evolution_AGN.param', 'r') as file:
    original_content = file.readlines()


TA = time.time()

N_lst = len(lst)

for i in range(N_lst):
  
  rkpc_i = lst[i][0]
  Lsh_i = lst[i][1]
  OutHDF5FileName = mainFunc(i, rkpc_i, Lsh_i)
  
  if i == 5:
    break


print(f'Elapsed time = {(time.time() - TA):.1f} sec')





