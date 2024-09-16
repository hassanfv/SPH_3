import h5py
import numpy as np



#------ Single Chimes Output -----
f = h5py.File(f'./SingleChimesRun.hdf5', 'r') #!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#TempEvol = f['TemperatureEvolution'][:]
AbundEvol = f['AbundanceEvolution'][:]     # (T, nH, Z, Elm, t)
#print(AbundEvol.shape)
Abchim = AbundEvol[0, 0, 0, :, 0]
print(Abchim)
#---------------------------------

#getInitIonState()


#!!!!! CHIMES - This is for 1/10 Solar !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! OVERRIDES ABOVE MASS FRACTION !!!!!!!!!!!!!!!
MassFrac = np.array([1.29e-03, 2.5306e-01, 2.07e-04, 8.36e-05, 5.49e-04, 1.41e-04, 5.91e-05, 6.83e-05, 4.09e-05, 6.44e-06, 1.1e-04]) # 0.1 Z_sun
print('MassFrac CHIMES = ', MassFrac)
print()
#----------------------

nHG = 10**np.arange(-1.0, 4.01, 0.2) #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
TempG = 10**np.arange(4.0, 7.21, 0.5)#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
N_nH = len(nHG)
N_T = len(TempG)
Npart = N_nH * N_T

nH_arr = np.zeros(Npart)
T_arr = np.zeros(Npart)
MassFrac_arr = np.zeros((Npart, 11))
coord_arr = np.zeros((Npart, 3))
initial_chemical_state = np.zeros((Npart, 157))

dist = 0.11 # kpc !!!!!!!!!!!!!!!!!!!!!!! TO BE VARIED !!!!!!!!

k = 0
for i in range(N_nH):
  for j in range(N_T):
    nH_arr[k] = nHG[i]
    T_arr[k] = TempG[j]
    MassFrac_arr[k, :] = MassFrac
    x = y = z = 1./np.sqrt(3.) * dist # we assume x = y = z for simplicity! dist = sqrt(x*x + y*y + z*z)
    coord_arr[k, :] = [x, y, z]
    initial_chemical_state[k, :] = Abchim
    k += 1


filename = 'hfvInput.hdf5'

with h5py.File(filename, 'w') as file:
    # Create a group named 'PartType0' in the HDF5 file to store datasets related to this particle type.
    part_type0 = file.create_group('PartType0')

    # Metallicity array Given as mass fractions relative to total in the order: 
    # All_metals, He, C, N, O, Ne, Mg, Si, S, Ca, Fe 
    part_type0.create_dataset('Metallicity', data=MassFrac_arr)
    
    nH = nHG
    part_type0.create_dataset('nHG_hfv', data=nH_arr)
    
    # Generate initial chemical states for each particle; it is for 157 species!
    part_type0.create_dataset('InitIonState_hfv', data=initial_chemical_state)
    
    # Temperature values for each particle.
    part_type0.create_dataset('TempG_hfv', data=T_arr)
    
    # 3D coordinates for each particle.
    part_type0.create_dataset('Coordinates', data=coord_arr)







