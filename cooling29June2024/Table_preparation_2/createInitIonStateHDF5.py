import h5py
import numpy as np

filename = 'sample_data.hdf5'

# Open a new HDF5 file to write data
with h5py.File(filename, 'w') as file:
    # Create a group named 'PartType0' in the HDF5 file to store datasets related to this particle type.
    part_type0 = file.create_group('PartType0')
    
    # Generate random data to simulate metallicities for N_part particles across 11 elements.
    metals = np.random.rand(N_part, 11)  # Random values for each metal for each particle
    part_type0.create_dataset('GFM_Metals', data=metals)
    
    # Generate random total metallicity values for each particle, non-logarithmic form.
    total_metallicity = np.random.rand(N_part) * 0.1  # Values between 0 and 0.1
    part_type0.create_dataset('GFM_Metallicity', data=total_metallicity)
    
    # Generate another set of random metallicity values for N_part particles across 11 elements.
    metallicity = np.random.rand(N_part, 11)  # Random values
    part_type0.create_dataset('Metallicity', data=metallicity)
    
    # Generate random hydrogen number density values for each particle.
    nH = np.random.rand(N_part) * 1e-4  # Values between 0 and 1e-4
    part_type0.create_dataset('nHG_hfv', data=nH)
    
    # Generate initial chemical states for each particle; using integer values as placeholders.
    initial_chemical_state = np.random.randint(0, 10, size=(N_part, 20))  # Integer states
    part_type0.create_dataset('InitIonState_hfv', data=initial_chemical_state)
    
    # Generate random temperature values for each particle.
    temperature = np.random.rand(N_part) * 1e4  # Temperatures between 0 and 10,000
    part_type0.create_dataset('TempG_hfv', data=temperature)
    
    # Generate random 3D coordinates for each particle.
    coordinates = np.random.rand(N_part, 3)  # 3D coordinates
    part_type0.create_dataset('Coordinates', data=coordinates)

