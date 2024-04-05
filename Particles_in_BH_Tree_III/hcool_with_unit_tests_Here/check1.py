
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Reference: Gnat & Sternberg - 2007 ----> https://cdsarc.cds.unistra.fr/viz-bin/cat/J/ApJS/168/213
df = pd.read_csv('table4a.dat', delim_whitespace=True, header=None)

T = df.iloc[:, 0].values
Lam = df.iloc[:, 3].values

print(df.head())


plt.plot(np.log10(T), np.log10(Lam), color = 'k')

plt.show()


