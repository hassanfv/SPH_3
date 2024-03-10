
import numpy as np
import pandas as pd


df = pd.read_csv('ngb.csv')
print(df)

counts = df.apply(lambda row: (row == -1).sum(), axis=1)

counts.to_csv('counts_of_negative_ones.csv', index = False, header = False)


