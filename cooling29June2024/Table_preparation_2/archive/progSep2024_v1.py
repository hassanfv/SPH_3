
import numpy as np
import pandas as pd
import pickle
import os
import time
import h5py
from mpi4py import MPI


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
LshG = np.arange(0.0, 2.51, 0.25) # it is in log10

lst = []
for rkpc in rkpcG:
  for Lsh in LshG:
    lst.append([rkpc, Lsh])
#--------------------------------------


#===== mainFunc
def mainFunc(nbeg, nend):

  for i in range(nbeg, nend):

    lst = inLists[i]

    OutFile = f'./grid_{i:03}_rkpc_{lst[0]:.2f}_Lsh_{lst[1]:.2f}.hdf5'

    Lsh = 10**lst[1] # Here the Lsh = min(Lsh, L_max) condition is implemented in "createInputList_Aug24.py"

    user_input = {
      "output_file": "   " + OutFile + '\n',
      "distance_to_AGN_kpc": f'{lst[0]:.3f}\n',
      "max_shield_length": " " + f'{Lsh:.3E}\n'}

    # Update the parameters in the file content
    updated_content = update_parameters(original_content, user_input)

    # Write the updated content to a new file
    updated_file_name = './grid_' + str(i) + '.param'
    with open(updated_file_name, 'w') as file:
        file.writelines(updated_content)

    #---- Executing CHIMES ----
    command = f"python3 chimes-driverx.py {updated_file_name}"
    os.system(command)
    #----

  return 1



# Read the initial param file.
with open('grid_noneq_evolution_AGN.param', 'r') as file:
    original_content = file.readlines()


with open('inputLists_Aug24.pkl', 'rb') as f:
  inLists = pickle.load(f)

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nCPUs = comm.Get_size()

N = len(inLists)

#------- used in MPI --------
count = N // nCPUs
remainder = N % nCPUs

if rank < remainder:
	nbeg = rank * (count + 1)
	nend = nbeg + count + 1
else:
	nbeg = rank * count + remainder
	nend = nbeg + count
#--------------------------

if rank == 0:
	T1 = time.time()
#--------------------------

local_res = mainFunc(nbeg, nend)

if rank == 0:
	tmp = local_res
	for i in range(1, nCPUs):
		tmp1 = comm.recv(source = i)
else:
	comm.send(local_res, dest = 0)
#----------------------------

if rank == 0:
  print(f'Elapsed time = ', time.time() - T1)




