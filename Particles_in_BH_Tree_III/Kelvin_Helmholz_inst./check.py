
import numpy as np


# Replace 'ngb_data.txt' with your actual filename
filename = 'ngb_data.txt'

# Read the file
with open(filename, 'r') as file:
    data = file.read()

# Convert the string data to a numpy array
# Assuming the data is separated by commas as in the C++ code above
ngb_array = np.fromstring(data, dtype=int, sep=',')

# Reshape the array if needed, assuming you know N and MAX_ngb
N = 1000000
MAX_ngb = 200
#ngb = ngb_array.reshape((N, MAX_ngb))
#print(np.sort(ngb))


contains_nan = np.isnan(ngb_array).any()

print("Array contains NaN:", contains_nan)


