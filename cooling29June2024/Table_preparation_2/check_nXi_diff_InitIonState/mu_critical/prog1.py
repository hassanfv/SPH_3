import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import h5py


#===== mainFunc
def mainFunc(nH_i, rkpc_i, Lsh_i, nYrs_i):

  OutFile = f'./rkpc_{rkpc_i:.2f}_Lsh_{np.log10(10**Lsh_i * pc_to_cm):.3f}.hdf5'
  OutFile_pkl = f'./rkpc_{rkpc_i:.2f}_Lsh_{np.log10(10**Lsh_i * pc_to_cm):.3f}.pkl'

  user_input = {
    "output_file": "   " + OutFile + '\n',
    "distance_to_AGN_kpc": f'{rkpc_i:.3f}\n',
    "max_shield_length": " " + f'{(10**Lsh_i * pc_to_cm):.3E}\n',
    "n_iterations": " " + f'{nYrs_i}',
    "log_nH_min": " " + f'{log_nH_i}',
    "log_nH_max": " " + f'{log_nH_i}'}

  # Update the parameters in the file content
  updated_content = update_parameters(original_content, user_input)

  # Write the updated content to a new file
  updated_file_name = f'rkpc_{rkpc_i:.2f}_Lsh_{np.log10(10**Lsh_i * pc_to_cm):.3f}.param'
  with open(updated_file_name, 'w') as file:
      file.writelines(updated_content)

  #---- Executing CHIMES ----
  #command = f"mpirun -np 1 python3 chimes-driver.py {updated_file_name}"
  command = f"mpirun -np 96 -x LD_LIBRARY_PATH=/path/to/install/dir/lib:$LD_LIBRARY_PATH python3 chimes-driver.py {updated_file_name}"
  print(command)
  os.system(command)
  #----
  
  #os.remove(updated_file_name)

  return OutFile, OutFile_pkl



#----------- Preparing the grid -------
rkpcG = np.arange(0.01, 0.52, 0.1)# it is in kpc
LshG = np.arange(0.0, 2.51, 0.25) # it is in log10 ---> We take Lsh in the range 1.0 pc up to ~300 pc.
#--------------------------------------

nYrsG = 










