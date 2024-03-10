
import numpy as np
import pandas as pd


df = pd.read_csv('counts_of_negative_ones.csv', header=None, index_col=False)
print(df)

val = df.values
val = val.flatten()
val = 200 - val

print(np.sort(val))

print(val[997859])

