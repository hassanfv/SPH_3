import h5py
import numpy as np



#------ Single Chimes Output -----
f = h5py.File(f'./SingleChimesRun.hdf5', 'r')
#TempEvol = f['TemperatureEvolution'][:]
AbundEvol = f['AbundanceEvolution'][:]     # (T, nH, Z, Elm, t)
#print(AbundEvol.shape)
Abchim = AbundEvol[0, 0, 0, :, 0]
print(Abchim)
#---------------------------------



'''
Amass = {
          "H": 1.008,    # Hydrogen
          "He": 4.0026,  # Helium
          "C": 12.011,   # Carbon
          "N": 14.007,   # Nitrogen
          "O": 15.999,   # Oxygen
          "Ne": 20.180,  # Neon
          "Mg": 24.305,  # Magnesium
          "Si": 28.085,  # Silicon
          "S": 32.06,    # Sulfur
          "Ca": 40.078,  # Calcium
          "Fe": 55.845   # Iron
}
'''


Amass = np.array([1.008, 4.0026, 12.011, 14.007, 15.999, 20.180, 24.305, 28.085, 32.06, 40.078, 55.845])
ElmAbund_solar = np.array([12.00, 10.93, 8.43, 7.83, 8.69, 7.93, 7.60, 7.51, 7.12, 6.34, 7.50])
ElmAbund_OneTenthSolar = np.array([12.00, 10.93, 7.43, 6.83, 7.69, 6.93, 6.60, 6.51, 6.12, 5.34, 6.50])

ElmMass = Amass * 10**(ElmAbund_solar - 12.0) # ARE YOU SURE WHICH ONE YOU CHOSE (solar or another)?????!!!!!!!!!!!!!!!!?????????
print(ElmMass)
print()

MassFrac = ElmMass / np.sum(ElmMass)
print('MassFrac Before = ', MassFrac)
print()
metal_mass = np.sum(ElmMass[2:]) / np.sum(ElmMass)
print('metal_mass = ', metal_mass)
print()
MassFrac[0] = metal_mass # as explained in snapshot_utils.py, the first value should be All metals mass fraction!
print('MassFrac After = ', MassFrac)
print()

#!!!!! CHIMES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! OVERRIDES ABOVE MASS FRACTION !!!!!!!!!!!!!!!
#MassFrac = np.array([0.0129, 0.2806, 2.07e-3, 8.36e-4, 5.49e-3, 1.41e-3, 5.91e-4, 6.83e-4, 4.09e-4, 6.44e-5, 1.1e-3]) # Original solar metallicity
MassFrac = np.array([1.29e-03, 2.5306e-01, 2.07e-04, 8.36e-05, 5.49e-04, 1.41e-04, 5.91e-05, 6.83e-05, 4.09e-05, 6.44e-06, 1.1e-04]) # 0.1 Z_sun
print('MassFrac CHIMES = ', MassFrac)
print()
#----------------------

nHG = 10**np.arange(-1.0, 4.01, 0.2)
TempG = 10**np.arange(4.0, 7.21, 0.5)
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

print(f'coord_arr = {coord_arr}')

print()
print(MassFrac_arr)

print('----------')
nt = np.where((nH_arr < 10000) & (nH_arr > 5000) & (T_arr < 6e5) & (T_arr > 1e5))
print(nt)
print()
print('nHHH = ', nH_arr[nt])
print()
print('TTT = ', T_arr[nt])
jj = 1
print(nH_arr[jj], T_arr[jj])




filename = 'hfvInput.hdf5'

# Open a new HDF5 file to write data
with h5py.File(filename, 'w') as file:
    # Create a group named 'PartType0' in the HDF5 file to store datasets related to this particle type.
    part_type0 = file.create_group('PartType0')

    # Generate set of metallicity values for N_part particles across 11 elements.
    # Metallicity array Given as mass fractions relative to total in the order: 
    # All_metals, He, C, N, O, Ne, Mg, Si, S, Ca, Fe 
    part_type0.create_dataset('Metallicity', data=MassFrac_arr)
    
    # Generate random hydrogen number density values for each particle.
    nH = nHG
    part_type0.create_dataset('nHG_hfv', data=nH_arr)
    
    # Generate initial chemical states for each particle; it is for 157 species!
    part_type0.create_dataset('InitIonState_hfv', data=initial_chemical_state)
    
    # Temperature values for each particle.
    part_type0.create_dataset('TempG_hfv', data=T_arr)
    
    # 3D coordinates for each particle.
    part_type0.create_dataset('Coordinates', data=coord_arr)







