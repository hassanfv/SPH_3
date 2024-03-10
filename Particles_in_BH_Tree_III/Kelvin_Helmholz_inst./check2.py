
import numpy as np
import pandas as pd


filename = 'ngb_data.txt'

with open(filename, 'r') as file:
    data = file.read()

ngb_array = np.fromstring(data, dtype=int, sep=',')

# Reshape the array if needed, assuming you know N and MAX_ngb
N = 1000000
MAX_ngb = 200
ngb = ngb_array.reshape((N, MAX_ngb))

df = pd.DataFrame(ngb)

df.to_csv('ngb.csv', index = False)

