import h5py
import numpy as np


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
ElmAbund = np.array([12.00, 10.93, 8.43, 7.83, 8.69, 7.93, 7.60, 7.51, 7.12, 6.34, 7.50])

ElmMass = Amass * 10**(ElmAbund - 12.0)
print(ElmMass)
print()

MassFrac = ElmMass / np.sum(ElmMass)
print('MassFrac = ', MassFrac)
print()

nHG = np.arange(-4.0, 4.01, 0.1)
TempG = np.arange(3.0, 10.21, 0.1)
N_nH = len(nHG)
N_T = len(TempG)
Npart = N_nH * N_T

nH_arr = np.zeros(Npart)
T_arr = np.zeros(Npart)
MassFrac_arr = np.zeros((Npart, 11))
coord_arr = np.zeros((Npart, 3))

dist = 0.11 # kpc !!!!!!!!!!!!!!!!!!!!!!! TO BE VARIED !!!!!!!!

k = 0
for i in range(N_nH):
  for j in range(N_T):
    nH_arr[k] = nHG[i]
    T_arr[k] = TempG[j]
    MassFrac_arr[k, :] = MassFrac
    r = dist**(1./3.)
    coord_arr[k, :] = [r, r, r]
    k += 1

print(coord_arr)

s()

filename = 'sample_data.hdf5'

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
    
    # Generate initial chemical states for each particle; it is for 157 species! Try to use CHIMES output for it ! THINK !!!!
    #initial_chemical_state = np.random.randint(0, 10, size=(N_part, 20))  # Integer states
    part_type0.create_dataset('InitIonState_hfv', data=initial_chemical_state)
    
    # Generate random temperature values for each particle.
    temperature = np.random.rand(N_part) * 1e4  # Temperatures between 0 and 10,000
    part_type0.create_dataset('TempG_hfv', data=T_arr)
    
    # Generate random 3D coordinates for each particle.
    coordinates = np.random.rand(N_part, 3)  # 3D coordinates
    part_type0.create_dataset('Coordinates', data=coord_arr)







