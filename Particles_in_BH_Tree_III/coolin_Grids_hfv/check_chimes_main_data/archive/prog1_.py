
import h5py
import matplotlib.pyplot as plt
import numpy as np

# Open the HDF5 file
with h5py.File('chimes_main_data.hdf5', 'r') as file:
    # Access the "cooling/rates" dataset
    if 'cooling/rates' in file:
        Rates = file['cooling/rates'][:]
        print("Shape of 'cooling/rates':", Rates.shape)
    else:
        print("Dataset 'cooling/rates' not found in the file.")
        Rates = None  # or handle the absence in another appropriate way

    # Access the "TableBins/Temperatures" dataset
    if 'TableBins/Temperatures' in file:
        Temp = file['TableBins/Temperatures'][:]
        print("Shape of 'TableBins/Temperatures':", Temp.shape)
    else:
        print("Dataset 'TableBins/Temperatures' not found in the file.")
        Temp = None  # or handle the absence in another appropriate way


plt.scatter(Temp, (Rates[14, :]))

plt.show()


